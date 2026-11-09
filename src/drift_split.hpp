#pragma once
struct Split {
 vector<int> ids;std::unique_ptr<Drift> head;std::unique_ptr<Bands> tail;
 Split(const vector<Fact>& f,int t,int g){
 vector<int> counts(g);for(auto e:f)++counts[e.key];
 for(int i=0;i<g;++i)if(counts[i])ids.push_back(i);
 std::sort(ids.begin(),ids.end(),[&](int a,int b){return counts[a]!=counts[b]?counts[a]>counts[b]:a<b;});
 size_t budget=(size_t(8)<<20)/(4*(size_t(t)+1))/32*32;
 ids.resize(std::min(ids.size(),std::min<size_t>(4096,budget)));ids.shrink_to_fit();std::sort(ids.begin(),ids.end());
 vector<int> local(g,-1);for(int i=0;i<(int)ids.size();++i)local[ids[i]]=i;
 vector<Fact> h,c;for(auto e:f){int i=local[e.key];if(i<0)c.push_back(e);else h.push_back({e.time,i});}
 head=std::make_unique<Drift>(h,t,int(ids.size()),true);tail=std::make_unique<Bands>(c,t,g,"band-raw");
 }
 i64 bytes(int mode=3)const{return sizeof(*this)+4LL*ids.capacity()+head->bytes(mode)+tail->bytes();}
 vector<Hit> query(Query q,int mode,bool share,DriftWork& dw,Metrics& w){
 auto h=head->query(q,mode,dw);for(auto& x:h)x.key=ids[x.key];
 if(share)return tail->seeded(q,w,h);
 auto c=tail->query(q,w);h.insert(h.end(),c.begin(),c.end());std::sort(h.begin(),h.end(),[](Hit x,Hit y){return std::abs(x.delta)!=std::abs(y.delta)?std::abs(x.delta)>std::abs(y.delta):x.key<y.key;});
 if(h.size()>size_t(q.k))h.resize(q.k);return h;}
};
