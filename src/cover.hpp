#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using std::vector;
using i64 = std::int64_t;
using Clock = std::chrono::steady_clock;
struct Fact {int time,key;};
struct Cell {int key,count;};
struct Query {int id,a,b,c,d,k;};
struct Hit {int key,delta;};
struct Metrics {int verified=0, passes=0, full=0, lists=0, ids=0, details=0;};
using Histogram=vector<Cell>;

inline void require(bool condition,const char* text){
    if(!condition) throw std::runtime_error(text);
}
inline int overlap(int lo,int hi,int begin,int end){
    return std::max(0,std::min(hi,end)-std::max(lo,begin));
}
struct Top {
    using Rank=std::pair<int,int>;
    std::priority_queue<Rank,vector<Rank>,std::greater<Rank>> heap;
    int k;
    explicit Top(int requested):k(requested){}
    void add(int key,int delta){
        if(delta==0) return;
        Rank rank{std::abs(delta),-key};
        if(int(heap.size())==k){
            if(rank<=heap.top()) return;
            heap.pop();
        }
        heap.push(rank);
    }
    int threshold() const {return int(heap.size())==k?heap.top().first:0;}
    template<class Score> vector<Hit> finish(Score score){
        vector<Hit> hits(heap.size());
        for(size_t i=hits.size();i-->0;){
            int key=-heap.top().second;heap.pop();hits[i]={key,score(key)};
        }
        return hits;
    }
};


struct History {
    vector<Cell> changes;
    vector<int> dense;
    int prefix(int end) const {
        if(!dense.empty()) return dense[end];
        auto it=std::lower_bound(changes.begin(),changes.end(),end,
                    [](Cell c,int t){return c.key<t;});
        return it==changes.begin()?0:(it-1)->count;
    }
    void finish(int bins){
        int total=0;
        for(auto& entry:changes){total+=entry.count;entry.count=total;}
        if(2*changes.size()>=size_t(bins+1)){
            dense.resize(bins+1);
            for(auto c:changes) dense[c.key+1]=c.count;
            for(int t=1;t<=bins;++t) dense[t]=std::max(dense[t],dense[t-1]);
            Histogram().swap(changes);
        }else changes.shrink_to_fit();
    }
};
struct RankedList {
    Histogram entries;
    void build(Histogram& values,int limit){
        auto earlier=[](Cell a,Cell b){return a.count!=b.count?a.count>b.count:a.key<b.key;};
        size_t n=std::min(values.size(),size_t(limit)+1);
        if(n<values.size()) std::nth_element(values.begin(),values.begin()+n,values.end(),earlier);
        entries.assign(values.begin(),values.begin()+n);
        std::sort(entries.begin(),entries.end(),earlier);
    }
    int tail(int n) const {return n<int(entries.size())?entries[n].count:0;}
};
struct NodeLists {RankedList raw,up,down;};
enum class Basis {Raw,Haar,Mixed};
#include "circuit.inc"

struct Index {
    int bins,groups,padded=1,cap,detail_cap;
    Basis basis;
    vector<NodeLists> tree;
    vector<History> histories;
    vector<int> seen;
    int epoch=0;
    Histogram frequency_order;
    bool bounded_completion=true;
    i64 history_bytes=0,summary_bytes=0;
    Index(const vector<Fact>& facts,int t,int g,int limit,Basis layout,int detail_limit=-1)
        :bins(t),groups(g),cap(limit),detail_cap(detail_limit<0?limit:detail_limit),basis(layout),histories(g),seen(g){
        while(padded<t) padded*=2;
        tree.resize(2*padded);
        vector<Histogram> leaves(padded);
        vector<int> counts(g),touched;
        size_t at=0;
        for(int time=0;time<t;++time){
            touched.clear();
            while(at<facts.size()&&facts[at].time==time){
                int key=facts[at++].key;
                if(counts[key]++==0) touched.push_back(key);
            }
            std::sort(touched.begin(),touched.end());
            for(int key:touched){
                leaves[time].push_back({key,counts[key]});
                histories[key].changes.push_back({time,counts[key]});
                counts[key]=0;
            }
        }
        require(at==facts.size(),"unsorted or invalid event time");
        for(auto& h:histories){
            h.finish(t);history_bytes+=sizeof(History)+h.dense.capacity()*sizeof(int)+h.changes.capacity()*sizeof(Cell);
        }
        frequency_order=build(1,leaves);
        std::sort(frequency_order.begin(),frequency_order.end(),[](Cell a,Cell b){
            return a.count!=b.count?a.count>b.count:a.key<b.key;
        });
        summary_bytes=tree.capacity()*sizeof(NodeLists);
        for(auto& n:tree) for(auto list:{&n.raw,&n.up,&n.down})
            summary_bytes+=list->entries.capacity()*sizeof(Cell);
    }
    Histogram build(int node,vector<Histogram>& leaves){
        Histogram counts;
        if(node>=padded) counts=std::move(leaves[node-padded]);
        else {
            auto left=build(2*node,leaves),right=build(2*node+1,leaves);
            Histogram up,down;size_t i=0,j=0;
            while(i<left.size()||j<right.size()){
                int key=std::min(i<left.size()?left[i].key:groups,j<right.size()?right[j].key:groups);
                int a=0,b=0;
                if(i<left.size()&&left[i].key==key) a=left[i++].count;
                if(j<right.size()&&right[j].key==key) b=right[j++].count;
                counts.push_back({key,a+b});
                if(basis!=Basis::Raw){
                    if(a>b) up.push_back({key,a-b});
                    if(b>a) down.push_back({key,b-a});
                }
            }
            tree[node].up.build(up,detail_cap);tree[node].down.build(down,detail_cap);
        }
        if(basis!=Basis::Haar||node==1){
            Histogram copy=counts;tree[node].raw.build(copy,cap);
        }
        return counts;
    }
    int score(int key,Query q) const {
        const auto& h=histories[key];
        return (h.prefix(q.d)-h.prefix(q.c))-(h.prefix(q.b)-h.prefix(q.a));
    }
    vector<Hit> query(Query q,Metrics& work){
        if(q.a==q.c&&q.b==q.d) return {};
        if(epoch==std::numeric_limits<int>::max()){std::fill(seen.begin(),seen.end(),0);epoch=0;}
        ++epoch;Top top(q.k);
        auto verify=[&](int key){
            if(seen[key]==epoch) return;
            seen[key]=epoch;++work.verified;top.add(key,score(key,q));
        };
        Circuit circuit(tree,padded,q);
        for(int depth=1;;depth=std::min(cap,4*depth)){
            ++work.passes;
            i64 worst=0;
            for(int sign:{1,-1}){
                worst=std::max(worst,circuit.plan(basis,sign,depth,detail_cap));
                auto visit=[&](const RankedList& list,bool detail){
                    if(list.entries.empty())return;
                    ++work.lists;work.details+=detail;
                    int n=std::min(detail?std::min(detail_cap,depth):depth,int(list.entries.size()));work.ids+=n;
                    for(int j=0;j<n;++j)verify(list.entries[j].key);
                };
                circuit.read(circuit.count-1,sign,false,visit);
            }
            if(worst==0 || worst/i64(padded)<top.threshold()) break;
            if(depth==cap){
                work.full=1;
                if(bounded_completion){
                    for(auto item:frequency_order){
                        if(item.count<top.threshold())break;
                        verify(item.key);
                    }
                } else for(int key=0;key<groups;++key) verify(key);
                break;
            }
        }
        return top.finish([&](int key){return score(key,q);});
    }
};
