#pragma once
#include "cover.hpp"
#include <cstdint>
#include <stdexcept>


struct DriftControlWork {
    i64 hot_scores=0, cold_scores=0, key_bounds=0, key_skips=0;
    i64 block_bounds=0, block_skips=0, lifetime_skips=0;
};

class DriftControls {
    struct Key { int id, freq; };
    struct Point { int t, prefix; };
    struct Term { int t, lo, hi, coefficient; };
    struct Entry { int id, score; };
    struct Better {
        bool operator()(Entry x, Entry y) const {
            int a=std::abs(x.score), b=std::abs(y.score);
            return a>b || (a==b && x.id<y.id);
        }
    };
    struct Best {
        int k;
        std::vector<Entry> queue;
        explicit Best(int count):k(count) { queue.reserve(count); }
        bool reject(i64 bound,int id) const {
            return static_cast<int>(queue.size())==k &&
                (bound<std::abs(queue.front().score) ||
                 (bound==std::abs(queue.front().score) && id>=queue.front().id));
        }
        void offer(int id,int value) {
            if (!value || reject(std::abs(value),id)) return;
            if(static_cast<int>(queue.size())==k) {
                std::pop_heap(queue.begin(),queue.end(),Better{});
                queue.pop_back();
            }
            queue.push_back({id,value});
            std::push_heap(queue.begin(),queue.end(),Better{});
        }
        std::vector<Hit> finish() {
            std::sort(queue.begin(),queue.end(),Better{});
            std::vector<Hit> result;
            result.reserve(queue.size());
            for(auto v:queue) result.push_back({v.id,v.score});
            return result;
        }
    };

    static constexpr int lanes=64;
    int T, mode_cap, hot=0, coarse=0, width=1, slabs=0;
    std::vector<Key> keys;
    std::vector<int> dense, offsets;
    std::vector<Point> history;
    std::vector<std::uint16_t> checkpoints;
    std::vector<unsigned char> shifts;
    std::vector<int> minimum,maximum,min_id;

    template<class V> static i64 mem(const V& v) {return i64(v.capacity())*sizeof(typename V::value_type);}
    int prefix(int rank,int t) const {
        if(t<=0) return 0;
        if(t>=T) return keys[rank].freq;
        int i=rank-hot;
        auto first=history.begin()+offsets[i], last=history.begin()+offsets[i+1];
        auto it=std::lower_bound(first,last,t,[](const Point& p,int time){return p.t<time;});
        return it==first ? 0 : (it-1)->prefix;
    }
    int score(int rank, const Query& q) const {
        int i=rank-hot, first=offsets[i], last=offsets[i+1];
        if(last-first<=8) {
            int previous=0, change=0;
            for(int p=first;p<last;++p) {
                auto h=history[p];
                int sign=int(q.c<=h.t && h.t<q.d)-int(q.a<=h.t && h.t<q.b);
                change+=sign*(h.prefix-previous);
                previous=h.prefix;
            }
            return change;
        }
        return (prefix(rank,q.d)-prefix(rank,q.c))-(prefix(rank,q.b)-prefix(rank,q.a));
    }
    int terms(const Query& q, std::array<Term,4>& out) const {
        const int ends[]={q.a,q.b,q.c,q.d}, signs[]={1,-1,-1,1};
        int count=0;
        for(int e=0;e<4;++e) {
            int t=ends[e];
            if(!t) continue;
            int i=0;
            while(i<count && out[i].t!=t) ++i;
            if(i==count) {
                int lo=t==T ? slabs : t/width, hi=t==T ? slabs : (t+width-1)/width;
                out[count++]={t,lo,hi,0};
            }
            out[i].coefficient+=signs[e];
        }
        return count;
    }
    template<class Bounds> static i64 magnitude(const std::array<Term,4>& ts,int n,Bounds get) {
        i64 lower=0,upper=0;
        for(int i=0;i<n;++i) {
            auto t=ts[i];
            int s=t.coefficient;
            if(!s) continue;
            auto range=get(t);
            lower+=i64(s)*(s>0?range.first:range.second);
            upper+=i64(s)*(s>0?range.second:range.first);
        }
        return std::max(upper,-lower);
    }
    void build_checkpoints() {
        width=std::max(1,int((i64(T)+31)/32));
        slabs=int((i64(T)+width-1)/width);
        while(hot+coarse<int(keys.size()) && keys[hot+coarse].freq>=std::max(1,(slabs+1)/2)) ++coarse;
        int blocks=(coarse+lanes-1)/lanes;
        shifts.resize(blocks);
        checkpoints.resize(size_t(slabs+1)*coarse);
        if(mode_cap>=2) {
            minimum.assign(size_t(slabs+1)*blocks,std::numeric_limits<int>::max());
            maximum.assign(size_t(slabs+1)*blocks,0);
            min_id.assign(blocks,std::numeric_limits<int>::max());
        }
        for(int i=0;i<coarse;++i) {
            int block=i/lanes, rank=i+hot;
            if(i%lanes==0) {
                while((keys[rank].freq>>shifts[block])>65535) ++shifts[block];
            }
            if(mode_cap>=2) min_id[block]=std::min(min_id[block],keys[rank].id);
            int p=offsets[i], end=offsets[i+1], cumulative=0;
            for(int j=0;j<=slabs;++j) {
                int time=static_cast<int>(std::min<i64>(i64(j)*width,T));
                while(p<end && history[p].t<time) cumulative=history[p++].prefix;
                checkpoints[size_t(j)*coarse+i]=static_cast<std::uint16_t>(cumulative>>shifts[block]);
                if(mode_cap>=2) {
                    size_t pos=size_t(j)*blocks+block;
                    minimum[pos]=std::min(minimum[pos],cumulative);
                    maximum[pos]=std::max(maximum[pos],cumulative);
                }
            }
        }
    }
public:


    DriftControls(const std::vector<Fact>& facts,int bins,int groups,int max_mode=2,int hot_limit=1024):
        T(bins),mode_cap(max_mode) {
        if(bins<0 || groups<0 || max_mode<0 || max_mode>2 || hot_limit<0)
            throw std::invalid_argument("DriftControls dimensions or mode");
        std::vector<int> freq(groups), occupied(groups), position(groups,-1);
        int previous=-1;
        for(auto f:facts) {
            if(f.time<previous || f.time<0 || f.time>=bins || f.key<0 || f.key>=groups)
                throw std::invalid_argument("DriftControls requires time-ordered valid facts");
            previous=f.time;
            if(freq[f.key]==std::numeric_limits<int>::max()) throw std::overflow_error("key frequency");
            ++freq[f.key];
            if(position[f.key]!=f.time) {++occupied[f.key]; position[f.key]=f.time;}
        }
        for(int id=0;id<groups;++id) if(freq[id]) keys.push_back({id,freq[id]});
        std::sort(keys.begin(),keys.end(),[](Key x,Key y){
            return x.freq!=y.freq?x.freq>y.freq:x.id<y.id;
        });
        int n=int(keys.size());
        hot=std::min(hot_limit,n);
        dense.assign((size_t(T)+1)*hot,0);
        offsets.resize(n-hot+1);
        for(int rank=0;rank<n;++rank) {
            int id=keys[rank].id;
            position[id]=rank;
            if(rank>=hot) {
                i64 end=i64(offsets[rank-hot])+occupied[id];
                if(end>std::numeric_limits<int>::max()) throw std::overflow_error("cold history offsets");
                offsets[rank-hot+1]=static_cast<int>(end);
            }
        }
        history.resize(offsets.back());
        std::vector<int> cursor(offsets.begin(),offsets.end());
        for(auto f:facts) {
            int rank=position[f.key];
            if(rank<hot) {++dense[(size_t(f.time)+1)*hot+rank];continue;}
            int i=rank-hot, &p=cursor[i];
            int cumul=p==offsets[i]?0:history[p-1].prefix;
            if(p>offsets[i] && history[p-1].t==f.time) ++history[p-1].prefix;
            else history[p++]={f.time,cumul+1};
        }
        for(int t=1;t<=T;++t) for(int j=0;j<hot;++j)
            dense[size_t(t)*hot+j]+=dense[size_t(t-1)*hot+j];
        if(mode_cap) build_checkpoints();
    }
    i64 bytes(int mode) const {
        if(mode<0 || mode>mode_cap) throw std::invalid_argument("unbuilt DriftControls mode");
        i64 total=mem(keys)+mem(dense)+mem(offsets)+mem(history);
        if(mode>=1) total+=mem(checkpoints)+mem(shifts);
        if(mode>=2) total+=mem(minimum)+mem(maximum)+mem(min_id);
        return total;
    }
    i64 bytes() const { return bytes(mode_cap); }
    std::vector<Hit> query(Query q,int mode,DriftControlWork& work) const {
        work={};
        if(mode<0 || mode>mode_cap) throw std::invalid_argument("unbuilt DriftControls mode");
        if(q.a<0 || q.a>q.b || q.b>T || q.c<0 || q.c>q.d || q.d>T)
            throw std::out_of_range("invalid intervals");
        int n=int(keys.size());
        if(q.k<=0 || !n || (q.a==q.c && q.b==q.d)) return {};
        Best best(std::min(q.k,n));
        if(hot) {
            auto ptr=[&](int t){return dense.data()+size_t(t)*hot;};
            const int *a=ptr(q.a), *b=ptr(q.b), *c=ptr(q.c), *d=ptr(q.d);
            for(int first=0;first<hot;first+=lanes) {
                if(best.reject(keys[first].freq,keys[first].id)) {
                    work.lifetime_skips=n-first; return best.finish();
                }
                int limit=std::min(hot-first,lanes);
                int scores[lanes];
                for(int j=0;j<limit;++j) {
                    int i=first+j;
                    scores[j]=(d[i]-c[i])-(b[i]-a[i]);
                }
                for(int j=0;j<limit;++j) best.offer(keys[first+j].id,scores[j]);
                work.hot_scores+=limit;
            }
        }
        std::array<Term,4> ts{};
        int count=mode?terms(q,ts):0, blocks=int(shifts.size());
        for(int i=0;i<n-hot;) {
            int rank=hot+i;
            if(best.reject(keys[rank].freq,keys[rank].id)) {
                work.lifetime_skips=n-rank;break;
            }
            int limit=mode && i<coarse ? std::min(coarse,((i/lanes)+1)*lanes) : n-hot;
            int block=i/lanes;
            if(mode==2 && i<coarse) {
                ++work.block_bounds;
                auto envelope=[&](Term t){
                    return std::make_pair(minimum[size_t(t.lo)*blocks+block],maximum[size_t(t.hi)*blocks+block]);
                };
                i64 bound=magnitude(ts,count,envelope);
                if(best.reject(bound,min_id[block])) {
                    work.block_skips+=limit-i;i=limit;continue;
                }
            }
            for(;i<limit;++i) {
                rank=hot+i;
                int id=keys[rank].id,f=keys[rank].freq;
                if(best.reject(f,id)) {
                    work.lifetime_skips=n-rank;return best.finish();
                }
                if(mode && i<coarse) {
                    ++work.key_bounds;
                    int shift=shifts[block];
                    auto envelope=[&](Term t){
                        if(t.t==T) return std::make_pair(i64(f),i64(f));
                        i64 low=i64(checkpoints[size_t(t.lo)*coarse+i])<<shift;
                        i64 high=i64(checkpoints[size_t(t.hi)*coarse+i])<<shift;
                        return std::make_pair(low,std::min<i64>(f,high+((i64(1)<<shift)-1)));
                    };
                    if(best.reject(magnitude(ts,count,envelope),id)) {++work.key_skips;continue;}
                }
                ++work.cold_scores;
                best.offer(id,score(rank,q));
            }
        }
        return best.finish();
    }
    std::vector<Hit> query(Query q,int mode) const {
        DriftControlWork work;
        return query(q,mode,work);
    }
};
