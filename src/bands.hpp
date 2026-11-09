#pragma once

struct Bands {
    struct Band {
        int maximum=0;
        vector<int> keys;
        std::unique_ptr<Dense> dense;
        std::unique_ptr<Index> index;
    };
    vector<Band> bands;
    vector<int> output_scores;
    bool direct_anchor,share_cutoff;
    i64 retained=0;
    Bands(const vector<Fact>& facts,int bins,int groups,const std::string& method)
        :output_scores(groups),direct_anchor(method!="band-noanchor"&&method!="band-raw"),
         share_cutoff(method!="band-noshare"){
        vector<int> frequency(groups),local(groups),bucket(groups);
        for(auto e:facts)++frequency[e.key];
        vector<Band> unsorted(16);vector<vector<Fact>> parts(16);
        for(int key=0;key<groups;++key)if(frequency[key]){
            int level=0;for(int f=frequency[key];f>=4;f/=4)++level;
            bucket[key]=level;auto& b=unsorted[level];
            local[key]=int(b.keys.size());b.keys.push_back(key);b.maximum=std::max(b.maximum,frequency[key]);
        }
        for(auto e:facts)parts[bucket[e.key]].push_back({e.time,local[e.key]});
        for(int level=15;level>=0;--level){
            auto& b=unsorted[level];if(b.keys.empty())continue;
            int n=int(b.keys.size());
            if(n<=256){
                b.dense=std::make_unique<Dense>(parts[level],bins,n);
                retained+=sizeof(Dense)+sizeof(int)*b.dense->table.capacity();
            }else{
                Basis basis=method=="band-raw"?Basis::Raw:Basis::Mixed;
                b.index=std::make_unique<Index>(parts[level],bins,n,basis==Basis::Raw?194:128,basis,32);
                const auto& x=*b.index;
                retained+=sizeof(Index)+x.summary_bytes+x.history_bytes+x.frequency_order.capacity()*sizeof(Cell)
                         +x.seen.capacity()*sizeof(int);
            }
            retained+=b.keys.capacity()*sizeof(int);
            bands.push_back(std::move(b));
            vector<Fact>().swap(parts[level]);
        }
        retained+=sizeof(*this)+bands.capacity()*sizeof(Band)+output_scores.capacity()*sizeof(int);
    }
    i64 bytes() const{return retained;}
    static int sibling(Query q,int padded){
        int w=q.b-q.a,lo=std::min(q.a,q.c);
        if(w<=0||w!=q.d-q.c||(w&(w-1))||std::abs(q.a-q.c)!=w||lo%(2*w))return 0;
        return (padded+lo)/(2*w);
    }
    vector<Hit> query(Query q,Metrics& work){return seeded(q,work,{});} vector<Hit> seeded(Query q,Metrics& work,const vector<Hit>& seed){
        if(q.a==q.c&&q.b==q.d)return {};
        Top global(q.k);
        auto emit=[&](int key,int score){output_scores[key]=score;global.add(key,score);};
        for(auto s:seed)emit(s.key,s.delta); for(auto& b:bands){
            int floor=share_cutoff?global.threshold():0;
            if(b.maximum<floor)break;
            if(b.dense){
                const auto& a=b.dense->table;size_t n=b.keys.size();
                for(size_t key=0;key<n;++key){
                    int score=(a[q.d*n+key]-a[q.c*n+key])-(a[q.b*n+key]-a[q.a*n+key]);
                    emit(b.keys[key],score);++work.verified;
                }
                continue;
            }
            Index& x=*b.index;
            int anchor=direct_anchor?sibling(q,x.padded):0;
            if(anchor){

                const auto& h=x.tree[anchor];int sign=q.a<q.c?-1:1;
                int bound=std::max(h.up.tail(x.detail_cap),h.down.tail(x.detail_cap));
                Top test(q.k);
                for(auto list:{&h.up,&h.down}){
                    int orientation=list==&h.up?sign:-sign;
                    int available=std::min(x.detail_cap,int(list->entries.size()));
                    for(int j=0;j<available;++j){
                        auto c=list->entries[j];test.add(c.key,c.count*orientation);
                    }
                    work.ids+=available;++work.lists;++work.details;
                }
                if(!bound || bound<std::max(floor,test.threshold())){
                    for(auto list:{&h.up,&h.down}){
                        int orientation=list==&h.up?sign:-sign;
                        int n=std::min(x.detail_cap,int(list->entries.size()));
                        for(int j=0;j<n;++j){auto c=list->entries[j];emit(b.keys[c.key],orientation*c.count);}
                    }
                    continue;
                }
            }
            if(x.epoch==std::numeric_limits<int>::max()){std::fill(x.seen.begin(),x.seen.end(),0);x.epoch=0;}
            ++x.epoch;Top local(q.k);
            auto cutoff=[&](){return share_cutoff?global.threshold():local.threshold();};
            auto verify=[&](int key){
                if(x.seen[key]==x.epoch)return;
                x.seen[key]=x.epoch;int score=x.score(key,q);++work.verified;
                emit(b.keys[key],score);if(!share_cutoff)local.add(key,score);
            };
            Circuit circuit(x.tree,x.padded,q);
            for(int depth=1;;depth=std::min(x.cap,depth*4)){
                ++work.passes;i64 bound=0;
                for(int sign:{1,-1}){
                    bound=std::max(bound,circuit.plan(x.basis,sign,depth,x.detail_cap));
                    auto read=[&](const RankedList& list,bool detail){
                        int limit=detail?std::min(x.detail_cap,depth):depth;
                        int n=std::min(limit,int(list.entries.size()));work.ids+=n;
                        if(n){++work.lists;work.details+=detail;}
                        for(int j=0;j<n;++j)verify(list.entries[j].key);
                    };
                    circuit.read(circuit.count-1,sign,false,read);
                }
                if(!bound||bound/x.padded<cutoff())break;
                if(depth==x.cap){
                    ++work.full;
                    for(auto item:x.frequency_order){if(item.count<cutoff())break;verify(item.key);}
                    break;
                }
            }
        }
        return global.finish([&](int key){return output_scores[key];});
    }
};
