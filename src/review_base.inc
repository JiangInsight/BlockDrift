#include "cover.hpp"
#include <memory>
#include <sstream>
#include <windows.h>
#include <psapi.h>

double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
struct Dense {
    int groups;vector<int> table;
    Dense(const vector<Fact>& facts,int bins,int g):groups(g),table(size_t(bins+1)*g){
        for(auto e:facts) ++table[size_t(e.time+1)*g+e.key];
        for(size_t i=g;i<table.size();++i) table[i]+=table[i-g];
    }
    vector<Hit> query(Query q){
        Top top(q.k);
        auto score=[&](int key){return (table[size_t(q.d)*groups+key]-table[size_t(q.c)*groups+key])-
                                      (table[size_t(q.b)*groups+key]-table[size_t(q.a)*groups+key]);};
        for(int key=0;key<groups;++key) top.add(key,score(key));
        return top.finish(score);
    }
};
struct Scan {
    const vector<Fact>& facts;vector<int> offset,count,seen,touched;int epoch=0;
    Scan(const vector<Fact>& f,int t,int g):facts(f),offset(t+1),count(g),seen(g){
        for(auto e:f) ++offset[e.time+1];
        std::partial_sum(offset.begin(),offset.end(),offset.begin());
    }
    vector<Hit> query(Query q){
        ++epoch;touched.clear();Top top(q.k);
        std::array<int,4> cut{q.a,q.b,q.c,q.d};std::sort(cut.begin(),cut.end());
        for(size_t j=1;j<cut.size();++j){
            int t=cut[j-1];
            int sign=int(q.c<=t&&t<q.d)-int(q.a<=t&&t<q.b);
            if(!sign) continue;
            for(int i=offset[t];i<offset[cut[j]];++i){
                int k=facts[i].key;
                if(seen[k]!=epoch){seen[k]=epoch;count[k]=0;touched.push_back(k);}
                count[k]+=sign;
            }
        }
        for(int k:touched) top.add(k,count[k]);
        return top.finish([&](int k){return count[k];});
    }
};


struct Lifetime {
    vector<History> histories;Histogram order;
    Lifetime(const vector<Fact>& facts,int t,int g):histories(g){
        vector<int> total(g);
        for(auto e:facts){
            ++total[e.key];auto& h=histories[e.key].changes;
            if(h.empty()||h.back().key!=e.time) h.push_back({e.time,1});else ++h.back().count;
        }
        for(int i=0;i<g;++i){histories[i].finish(t);order.push_back({i,total[i]});}
        std::sort(order.begin(),order.end(),[](Cell a,Cell b){return a.count!=b.count?a.count>b.count:a.key<b.key;});
    }
    vector<Hit> query(Query q,Metrics& work){
        Top top(q.k);
        auto score=[&](int k){const auto& h=histories[k];return (h.prefix(q.d)-h.prefix(q.c))-(h.prefix(q.b)-h.prefix(q.a));};
        for(auto item:order){
            if(item.count<top.threshold()) break;
            ++work.verified;top.add(item.key,score(item.key));
        }
        return top.finish(score);
    }
};
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#endif

using Clock = std::chrono::steady_clock;
using Count = std::int64_t;
using Pos = std::uint32_t;

static int popcount(std::uint64_t value) {
#ifdef _MSC_VER
    return static_cast<int>(__popcnt64(value));
#else
    return __builtin_popcountll(value);
#endif
}

struct Level {
    std::vector<std::uint64_t> bits;
    std::vector<Pos> prefix;
    Pos zeros;
    Pos ones_before(Pos p) const {
        const auto word = p >> 6;
        const auto shift = p & 63;
        const auto mask = (std::uint64_t{1} << shift) - 1;
        return prefix[word] + (shift ? popcount(bits[word] & mask) : 0);
    }
};

struct Interval {
    Pos first{}, end{};
    int sign{};
};
struct Node {
    std::array<Interval, 2> portions{};
    Count bound{};
    std::uint32_t min_key{};
    int depth{};
};
struct NextNode {
    bool operator()(const Node& x, const Node& y) const {
        return x.bound != y.bound ? x.bound < y.bound : x.min_key > y.min_key;
    }
};
struct Answer {
    std::uint32_t key;
    Count delta;
};
struct Better {
    bool operator()(const Answer& x, const Answer& y) const {
        return std::abs(x.delta) != std::abs(y.delta) ? std::abs(x.delta) > std::abs(y.delta) : x.key < y.key;
    }
};
struct Counters { std::uint64_t popped{}, expanded{}, leaves{}; };

class Wavelet {
    std::vector<Level> levels;
    std::vector<Pos> offsets;
    int width = 0;

    static Count bound(const Node& s) {
        Count positive = 0, negative = 0;
        for (const auto& r : s.portions) {
            if (r.sign > 0) positive += r.end - r.first;
            if (r.sign < 0) negative += r.end - r.first;
        }
        return std::max(positive, negative);
    }
    static bool can_improve(const Node& n, std::size_t k,
            const std::priority_queue<Answer,std::vector<Answer>,Better>& best) {
        if (!n.bound) return false;
        if (best.size() < k) return true;
        const auto threshold = std::abs(best.top().delta);
        return n.bound > threshold || (n.bound == threshold && n.min_key < best.top().key);
    }
public:
    Count packed_read_ns = 0;
    Wavelet(const std::string& path, int bins, int groups) : offsets(static_cast<std::size_t>(bins) + 1) {
        if (bins < 0 || groups < 0) throw std::runtime_error("negative dimensions");
        auto read_start = Clock::now();
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        if (!in) throw std::runtime_error("cannot open packed stream");
        const auto bytes = static_cast<std::streamoff>(in.tellg());
        if (bytes < 0 || bytes % 8 || std::uint64_t(bytes/8) > std::numeric_limits<Pos>::max())
            throw std::runtime_error("invalid packed stream or 32-bit row limit exceeded");
        in.seekg(0);
        std::vector<Pos> keys;
        keys.reserve(static_cast<std::size_t>(bytes/8));
        int last = -1;
        std::array<std::int32_t, 16384> chunk;
        while (in.read(reinterpret_cast<char*>(chunk.data()), sizeof(chunk)) || in.gcount() > 0) {
            auto count = static_cast<std::size_t>(in.gcount()) / sizeof(chunk[0]);
            for (std::size_t i = 0; i < count; i += 2) {
                int t = chunk[i], v = chunk[i+1];
                if (t < last || t < 0 || t >= bins || v < 0 || v >= groups)
                    throw std::runtime_error("out-of-domain or unsorted input");
                ++offsets[t+1];
                keys.push_back(static_cast<Pos>(v));
                last = t;
            }
        }
        packed_read_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-read_start).count();
        for (std::size_t i = 1; i < offsets.size(); ++i) offsets[i] += offsets[i-1];
        for (auto max = groups ? static_cast<Pos>(groups-1) : 0; max; max >>= 1) ++width;
        std::vector<Pos> next(keys.size());
        for (int shift = width-1; shift >= 0; --shift) {
            Level l;
            const auto words = (keys.size()+63)/64;
            l.bits.resize(words);
            l.prefix.resize(words+1);
            for (std::size_t i = 0; i < keys.size(); ++i)
                if ((keys[i] >> shift) & 1u) l.bits[i/64] |= std::uint64_t{1} << (i%64);
            for (std::size_t i = 0; i < words; ++i) l.prefix[i+1] = l.prefix[i] + popcount(l.bits[i]);
            l.zeros = static_cast<Pos>(keys.size()) - l.prefix.back();
            std::size_t zero = 0, one = l.zeros;
            for (const auto key : keys) {
                if ((key >> shift) & 1u) next[one++] = key;
                else next[zero++] = key;
            }
            levels.push_back(std::move(l));
            keys.swap(next);
        }
    }
    std::size_t bytes() const {
        std::size_t s = sizeof(*this) + offsets.capacity()*sizeof(Pos) + levels.capacity()*sizeof(Level);
        for (const auto& l : levels) s += l.bits.capacity()*sizeof(std::uint64_t) + l.prefix.capacity()*sizeof(Pos);
        return s;
    }
    std::vector<Answer> query(int a, int b, int c, int d, int k, Counters& counters) const {
        if (a<0 || b<a || c<0 || d<c || std::max(b,d)>=static_cast<int>(offsets.size()) || k<0)
            throw std::runtime_error("out-of-domain query");
        if (!k) return {};
        std::array<int,4> v{a,b,c,d};
        std::sort(v.begin(),v.end());
        Node root;
        int used=0;
        for (int i=0;i<3;++i) {
            if (v[i]==v[i+1]) continue;
            int sign=(v[i]>=c && v[i]<d)-(v[i]>=a && v[i]<b);
            if (sign) root.portions.at(used++)={offsets[v[i]],offsets[v[i+1]],sign};
        }
        root.bound=bound(root);
        std::priority_queue<Node,std::vector<Node>,NextNode> open;
        std::priority_queue<Answer,std::vector<Answer>,Better> best;
        if (root.bound) open.push(root);
        while (!open.empty()) {
            const auto n=open.top();
            if (!can_improve(n,static_cast<std::size_t>(k),best)) break;
            open.pop();
            ++counters.popped;
            if (n.depth==width) {
                ++counters.leaves;
                Count delta=0;
                for (const auto& r:n.portions) delta += (Count(r.end)-r.first)*r.sign;
                if (!delta) continue;
                const Answer row{n.min_key,delta};
                if (best.size()<static_cast<std::size_t>(k)) best.push(row);
                else if (Better{}(row,best.top())) { best.pop(); best.push(row); }
                continue;
            }
            ++counters.expanded;
            const auto& l=levels[n.depth];
            std::array<Node,2> children{n,n};
            for (int j=0;j<2;++j) {
                children[j].depth=n.depth+1;
                if (j) children[j].min_key |= Pos{1} << (width-n.depth-1);
            }
            for (int i=0;i<used;++i) {
                const auto r=n.portions[i];
                auto l1=l.ones_before(r.first), h1=l.ones_before(r.end);
                children[0].portions[i]={r.first-l1,r.end-h1,r.sign};
                children[1].portions[i]={l.zeros+l1,l.zeros+h1,r.sign};
            }
            for (auto& child:children) {
                child.bound=bound(child);
                if (can_improve(child,static_cast<std::size_t>(k),best)) open.push(child);
            }
        }
        std::vector<Answer> result;
        while (!best.empty()) {result.push_back(best.top()); best.pop();}
        std::reverse(result.begin(), result.end());
        return result;
    }
};


#include "bands.hpp"
#include <map>

struct Materialization {
    i64 build_ns, retained_bytes;
};

template<class Factory, class Size>
auto materialize(std::map<std::string,Materialization>& records, const std::string& name,
                 Factory make, Size bytes) {
    auto begin=Clock::now();
    auto index=make();
    auto stop=Clock::now();
    records[name]={std::chrono::duration_cast<std::chrono::nanoseconds>(stop-begin).count(),i64(bytes(*index))};
    return index;
}

static std::map<int,vector<Hit>> read_oracle(const char* path) {
    std::ifstream file(path);
    require(bool(file),"missing oracle");
    std::map<int,vector<Hit>> expected;
    std::string line;
    std::getline(file,line);
    while(std::getline(file,line)) {
        std::replace(line.begin(),line.end(),',',' ');
        std::istringstream row(line);int id,rank;Hit hit;
        require(bool(row>>id>>rank>>hit.key>>hit.delta),"malformed oracle");
        auto& prior=expected[id];
        require(rank==int(prior.size())+1 && hit.delta!=0,"incorrect oracle rank or zero");
        prior.push_back(hit);
    }
    return expected;
}

template<class Function>
i64 timed_check(Function run, const vector<Hit>& expected, int qid, const std::string& method) {
    auto start=Clock::now();
    const auto result=run();
    auto stop=Clock::now();
    auto fail=[&]{throw std::runtime_error("answer mismatch for "+method+", qid "+std::to_string(qid));};
    if(result.size()!=expected.size()) fail();
    for(size_t i=0;i<result.size();++i)
        if(i64(result[i].key)!=expected[i].key || result[i].delta!=expected[i].delta) fail();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count();
}


class LayerMatrix {
    vector<Level> levels;
    vector<Pos> offsets;
    int width=0;
    static Count bound(const Node& n){
        Count positive=0,negative=0;
        for(auto& r:n.portions){
            if(r.sign>0)positive+=r.end-r.first;
            if(r.sign<0)negative+=r.end-r.first;
        }
        return std::max(positive,negative);
    }
public:
    LayerMatrix(const vector<Fact>& facts,int bins,int groups):offsets(bins+1) {
        vector<Pos> keys;keys.reserve(facts.size());
        for(auto& f:facts){++offsets[f.time+1];keys.push_back(f.key);}
        std::partial_sum(offsets.begin(),offsets.end(),offsets.begin());
        for(auto x=groups-1;x;x>>=1)++width;
        vector<Pos> next(keys.size());
        for(int shift=width-1;shift>=0;--shift){
            Level l;const auto words=(keys.size()+63)/64;
            l.bits.resize(words);l.prefix.resize(words+1);
            for(size_t i=0;i<keys.size();++i)if((keys[i]>>shift)&1u)l.bits[i/64]|=std::uint64_t{1}<<(i%64);
            for(size_t i=0;i<words;++i)l.prefix[i+1]=l.prefix[i]+popcount(l.bits[i]);
            l.zeros=Pos(keys.size())-l.prefix.back();
            size_t zero=0,one=l.zeros;
            for(Pos key:keys){if((key>>shift)&1u)next[one++]=key;else next[zero++]=key;}
            levels.push_back(std::move(l));keys.swap(next);
        }
    }
    size_t bytes()const{
        size_t n=sizeof(*this)+offsets.capacity()*sizeof(Pos)+levels.capacity()*sizeof(Level);
        for(const auto& l:levels)n+=l.bits.capacity()*sizeof(std::uint64_t)+l.prefix.capacity()*sizeof(Pos);
        return n;
    }
    template<class Emit>
    void search(Query q,const vector<int>& ids,Top& global,Emit emit,Counters& counter)const{
        std::array<int,4> cut{q.a,q.b,q.c,q.d};std::sort(cut.begin(),cut.end());
        Node root;int used=0;
        for(int i=0;i<3;++i){
            int t=cut[i];if(t==cut[i+1])continue;
            int sign=(t>=q.c&&t<q.d)-(t>=q.a&&t<q.b);
            if(sign)root.portions.at(used++)={offsets[t],offsets[cut[i+1]],sign};
        }
        auto improves=[&](const Node& n){
            if(!n.bound)return false;
            if(int(global.heap.size())<global.k)return true;
            const auto& worst=global.heap.top();
            return n.bound>worst.first || (n.bound==worst.first && ids[n.min_key]<-worst.second);
        };
        root.bound=bound(root);
        std::priority_queue<Node,vector<Node>,NextNode> open;
        if(improves(root))open.push(root);
        while(!open.empty()){
            const Node n=open.top();
            if(!improves(n))break;
            open.pop();++counter.popped;
            if(n.depth==width){
                ++counter.leaves;
                int delta=0;
                for(const auto& r:n.portions)delta+=(int(r.end)-int(r.first))*r.sign;
                if(delta)emit(ids[n.min_key],delta);
                continue;
            }
            ++counter.expanded;
            const auto& level=levels[n.depth];
            std::array<Node,2> kids{n,n};
            kids[0].depth=kids[1].depth=n.depth+1;
            kids[1].min_key|=Pos{1}<<(width-1-n.depth);
            for(int j=0;j<used;++j){
                auto r=n.portions[j];auto a=level.ones_before(r.first),b=level.ones_before(r.end);
                kids[0].portions[j]={r.first-a,r.end-b,r.sign};
                kids[1].portions[j]={level.zeros+a,level.zeros+b,r.sign};
            }
            for(auto& kid:kids){kid.bound=bound(kid);if(improves(kid))open.push(kid);}
        }
    }
};

struct Transfer {
    struct Band {int maximum=0;vector<int> keys;std::unique_ptr<Dense> dense;std::unique_ptr<LayerMatrix> matrix;std::unique_ptr<Lifetime> lifetime;};
    vector<Band> bands;
    vector<int> output;
    size_t retained;
    Transfer(const vector<Fact>& facts,int bins,int groups,bool matrix):output(groups),retained(0){
        vector<int> frequency(groups),local(groups),bucket(groups);
        for(auto e:facts)++frequency[e.key];
        vector<Band> all(16);vector<vector<Fact>> parts(16);
        for(int k=0;k<groups;++k)if(frequency[k]){
            int level=0;for(int f=frequency[k];f>=4;f/=4)++level;
            bucket[k]=level;auto& b=all[level];b.maximum=std::max(b.maximum,frequency[k]);
            local[k]=int(b.keys.size());b.keys.push_back(k);
        }
        for(auto e:facts)parts[bucket[e.key]].push_back({e.time,local[e.key]});
        for(int i=15;i>=0;--i){
            auto& b=all[i];if(b.keys.empty())continue;
            int n=int(b.keys.size());
            if(n<=256){
                b.dense=std::make_unique<Dense>(parts[i],bins,n);
                retained+=sizeof(Dense)+sizeof(int)*b.dense->table.capacity();
            }else if(matrix){
                b.matrix=std::make_unique<LayerMatrix>(parts[i],bins,n);
                retained+=b.matrix->bytes();
            }else{
                b.lifetime=std::make_unique<Lifetime>(parts[i],bins,n);
                const auto& x=*b.lifetime;
                retained+=sizeof(Lifetime)+x.histories.capacity()*sizeof(History)+x.order.capacity()*sizeof(Cell);
                for(auto& h:x.histories)retained+=h.dense.capacity()*sizeof(int)+h.changes.capacity()*sizeof(Cell);
            }
            retained+=b.keys.capacity()*sizeof(int);
            bands.push_back(std::move(b));vector<Fact>().swap(parts[i]);
        }
        retained+=sizeof(*this)+bands.capacity()*sizeof(Band)+output.capacity()*sizeof(int);
    }
    vector<Hit> query(Query q,Metrics& m,Counters& c){
        if(q.a==q.c&&q.b==q.d)return {};
        Top global(q.k);
        auto emit=[&](int k,int d){output[k]=d;global.add(k,d);};
        for(auto& b:bands){
            if(b.maximum<global.threshold())break;
            if(b.dense){
                const auto& a=b.dense->table;size_t n=b.keys.size();
                for(size_t k=0;k<n;++k){
                    auto d=(a[q.d*n+k]-a[q.c*n+k])-(a[q.b*n+k]-a[q.a*n+k]);
                    emit(b.keys[k],d);++m.verified;
                }
            }else if(b.matrix)b.matrix->search(q,b.keys,global,emit,c);
            else {
                const auto& x=*b.lifetime;
                for(const auto& item:x.order){
                    int floor=global.threshold(),gid=b.keys[item.key];
                    if(item.count<floor)break;
                    if(floor>0 && item.count==floor && gid>=-global.heap.top().second)break;
                    const auto& h=x.histories[item.key];
                    const int d=(h.prefix(q.d)-h.prefix(q.c))-(h.prefix(q.b)-h.prefix(q.a));
                    emit(gid,d);++m.verified;
                }
            }
        }
        return global.finish([&](int k){return output[k];});
    }
};
