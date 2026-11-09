#pragma once

struct DriftWork {int seen=0,box=0,pair=0,scored=0,cold=0; i64 cells=0;};
struct Drift {
 static constexpr int B=32,K=16;
 int T,G,H; vector<int> order,rank,freq,offset,pos,dense,mu,lo,hi,edge,out,cell_offset,delta,stamp; vector<Cell> cold_cells;int epoch=0;
 Drift(const vector<Fact>& f,int t,int g,bool all_hot=false):T(t),G(g),rank(g),freq(g),out(g){
  for(auto v:f)++freq[v.key];
  for(int i=0;i<G;++i)if(freq[i])order.push_back(i);
  std::sort(order.begin(),order.end(),[&](int x,int y){return freq[x]!=freq[y]?freq[x]>freq[y]:x<y;});
  int n=int(order.size());for(int i=0;i<n;++i)rank[order[i]]=i;
  H=std::min(n,int(std::min<size_t>(32768,std::max<size_t>(256,2*f.size()/(T+1))))/B*B);
  if(all_hot)H=n; int nb=(H+B-1)/B;
  dense.assign(size_t(nb)*B*(T+1),0);mu.resize(size_t(nb)*(T+1));lo.resize(nb*K);hi.resize(nb*K);edge.assign(size_t(nb)*K*K,INT_MIN);
  offset.assign(n+1,0);for(int i=H;i<n;++i)offset[i+1]=offset[i]+freq[order[i]];
  pos.resize(offset.back());auto wr=offset;
  for(auto v:f){int r=rank[v.key];if(r<H)++dense[(size_t(r/B)*(T+1)+v.time+1)*B+r%B];else pos[wr[r]++]=v.time;}
  cell_offset.resize(T+1);delta.resize(n-H);stamp.resize(n-H,-1);vector<int> touched;
  size_t fp=0;for(int t0=0;t0<T;++t0){touched.clear();while(fp<f.size()&&f[fp].time==t0){int r=rank[f[fp++].key]-H;if(r>=0){if(!delta[r]++)touched.push_back(r);}}for(int r:touched){cold_cells.push_back({r,delta[r]});delta[r]=0;}cell_offset[t0+1]=int(cold_cells.size());}
  cold_cells.shrink_to_fit();
  for(int b=0;b<nb;++b){
   int len=std::min(B,H-b*B);auto a=&dense[size_t(b)*(T+1)*B];auto m=&mu[size_t(b)*(T+1)];
   int low[B][K],high[B][K];
   for(int t0=0;t0<K;++t0){lo[b*K+t0]=INT_MAX;hi[b*K+t0]=INT_MIN;for(int j=0;j<len;++j){low[j][t0]=INT_MAX;high[j][t0]=INT_MIN;}}
   for(int ti=0;ti<=T;++ti){
    auto row=a+size_t(ti)*B;i64 sum=0;
    for(int j=0;j<len;++j){if(ti)row[j]+=row[j-B];sum+=row[j];}
    m[ti]=int(sum/len);int tile=tag(ti);
    for(int j=0;j<len;++j){int z=row[j]-m[ti];low[j][tile]=std::min(low[j][tile],z);high[j][tile]=std::max(high[j][tile],z);
      lo[b*K+tile]=std::min(lo[b*K+tile],z);hi[b*K+tile]=std::max(hi[b*K+tile],z);}
   }
   for(int i=0;i<K;++i)for(int k=0;k<K;++k)if(low[0][i]!=INT_MAX&&low[0][k]!=INT_MAX)
    for(int j=0;j<len;++j)edge[(b*K+i)*K+k]=std::max(edge[(b*K+i)*K+k],high[j][i]-low[j][k]);
  }
 }
 int tag(int t)const{return int(i64(t)*K/(T+1));}
 i64 bytes(int mode=3)const{return sizeof(*this)+cold_cells.capacity()*sizeof(Cell)+4LL*(cell_offset.capacity()+delta.capacity()+stamp.capacity())+4LL*(order.capacity()+rank.capacity()+freq.capacity()+offset.capacity()+pos.capacity()+dense.capacity()+out.capacity()+(mode?mu.capacity()+lo.capacity()+hi.capacity():0)+(mode>=2?edge.capacity():0));}
 vector<Hit> query(Query q,int mode,DriftWork& w){
  if(q.k<=0||(q.a==q.c&&q.b==q.d))return {};
  int p[2]={q.a,q.d},m[2]={q.b,q.c};int nr=2;
  for(int i=0;i<nr;){bool same=false;for(int j=0;j<nr;++j)if(p[i]==m[j]){p[i]=p[nr-1];m[j]=m[nr-1];--nr;same=true;break;}if(!same)++i;}
  if(!nr)return {};Top top(q.k);
  auto offer=[&](int key,int s){if(s && std::abs(s)>=top.threshold()){out[key]=s;top.add(key,s);}};
  const int ap=tag(p[0]),am=tag(m[0]),bp=nr==2?tag(p[1]):0,bm=nr==2?tag(m[1]):0;
  for(int i=0;i<H;i+=B){
   int cut=top.threshold();if(freq[order[i]]<cut)break;int b=i/B;int n=std::min(B,H-i);
   if(i&&mode){++w.seen;auto mean=&mu[size_t(b)*(T+1)];i64 base=i64(mean[p[0]])-mean[m[0]];if(nr==2)base+=i64(mean[p[1]])-mean[m[1]];
    auto l=&lo[b*K],h=&hi[b*K];i64 up=base+h[ap]-l[am],dn=-base+h[am]-l[ap];if(nr==2){up+=i64(h[bp])-l[bm];dn+=i64(h[bm])-l[bp];}
    auto dead=[&](i64 u,i64 v){i64 z=std::max(u,v);return z==0||z<cut;};if(dead(up,dn)){++w.box;continue;}
    if(mode>=2){auto e=&edge[b*K*K];i64 u=e[ap*K+am],v=e[am*K+ap];if(nr==2){u+=e[bp*K+bm];v+=e[bm*K+bp];if(mode==3){u=std::min(u,i64(e[ap*K+bm])+e[bp*K+am]);v=std::min(v,i64(e[am*K+bp])+e[bm*K+ap]);}}
      if(dead(base+u,-base+v)){++w.pair;continue;}}
   }
   auto a=&dense[size_t(b)*(T+1)*B];auto x=a+q.a*B,y=a+q.b*B,z=a+q.c*B,t=a+q.d*B;
   int scores[B];for(int j=0;j<n;++j)scores[j]=(t[j]-z[j])-(y[j]-x[j]);
   w.scored+=n;for(int j=0;j<n;++j)if(std::abs(scores[j])>=cut)offer(order[i+j],scores[j]);
  }
  int cap=top.threshold(),endrank=int(order.size());if(cap&&H<endrank){int first=H,hi=endrank,last=-top.heap.top().second;while(first<hi){int mid=first+(hi-first)/2,k=order[mid];if(freq[k]>cap||(freq[k]==cap&&k<last))first=mid+1;else hi=mid;}endrank=first;}
  if(H<endrank){
   int ts[4]={q.a,q.b,q.c,q.d};std::sort(ts,ts+4);
   int sign[3];i64 nc=0;for(int j=0;j<3;++j){int t=ts[j];sign[j]=int(t>=q.c&&t<q.d)-int(t>=q.a&&t<q.b);if(sign[j])nc+=cell_offset[ts[j+1]]-cell_offset[t];}
   if(nc<2LL*(endrank-H)){
    if(epoch==INT_MAX){std::fill(stamp.begin(),stamp.end(),-1);epoch=0;}++epoch;vector<int> touched;
    for(int j=0;j<3;++j)if(sign[j])for(int z=cell_offset[ts[j]];z<cell_offset[ts[j+1]];++z){auto v=cold_cells[z];++w.cells;if(v.key>=endrank-H)continue;
     if(stamp[v.key]!=epoch){stamp[v.key]=epoch;delta[v.key]=0;touched.push_back(v.key);}delta[v.key]+=sign[j]*v.count;}
    for(int r:touched){offer(order[H+r],delta[r]);++w.cold;}
   }else for(int i=H;i<endrank;++i){int key=order[i],cut=top.threshold();if(freq[key]<cut||(cut>0&&freq[key]==cut&&key>=-top.heap.top().second))break;
    auto start=pos.data()+offset[i],end=pos.data()+offset[i+1];auto F=[&](int t){return std::lower_bound(start,end,t);};
    int v=int((F(q.d)-F(q.c))-(F(q.b)-F(q.a)));++w.cold;offer(key,v);
   }
  }
  return top.finish([&](int key){return out[key];});
 }
};
