struct Circuit {
    struct CellPlan {
        const NodeLists* lists; int left,right;
        i64 mean,contrast;
        i64 detail_cost=0,best=0;
        bool choose_basis=false;
    };
    const vector<NodeLists>& tree;
    const int padded;
    Query q;
    std::array<CellPlan,192> cells;
    int count=0;
    Circuit(const vector<NodeLists>& t,int p,Query query):tree(t),padded(p),q(query){fill(1,0,p);}
    i64 weight(int lo,int hi) const {
        return overlap(lo,hi,q.c,q.d)-overlap(lo,hi,q.a,q.b);
    }
    void fill(int node,int lo,int hi){
        int width=hi-lo; i64 mass=weight(lo,hi),contrast=0;
        int left=-1,right=-1;
        bool cut=false;
        for(int endpoint:{q.a,q.b,q.c,q.d})cut|=lo<endpoint&&endpoint<hi;
        if(cut){
            int mid=lo+width/2;
            fill(2*node,lo,mid);left=count-1;
            fill(2*node+1,mid,hi);right=count-1;
            contrast=(2*weight(lo,mid)-mass)*(padded/width);
        }
        require(count<int(cells.size()),"query circuit bound");
        cells[count++]={&tree[node],left,right,mass*(padded/width),contrast};
    }
    i64 plan(Basis basis,int sign,int depth,int detail_limit=std::numeric_limits<int>::max()){
        for(int i=0;i<count;++i){
            auto& p=cells[i]; const auto& n=*p.lists;
            i64 mean=sign*p.mean;
            p.best=mean>0?mean*n.raw.tail(depth):0;
            p.detail_cost=0;p.choose_basis=false;
            if(p.left<0)continue;
            auto& a=cells[p.left];auto& b=cells[p.right];i64 split=a.best+b.best;
            if(basis==Basis::Raw){p.best=split;continue;}
            i64 contrast=sign*p.contrast;
            p.detail_cost=a.detail_cost+b.detail_cost+
                std::abs(contrast)*(contrast>0?n.up:n.down).tail(std::min(depth,detail_limit));
            p.best+=p.detail_cost;
            p.choose_basis=(basis==Basis::Haar||p.best<split);
            if(!p.choose_basis)p.best=split;
        }
        return cells[count-1].best;
    }
    template<class Reader>
    void read(int ordinal,int sign,bool details_only,Reader& visit) const {
        const auto& p=cells[ordinal];
        if(details_only||p.choose_basis||p.left<0){
            if(!details_only&&sign*p.mean>0)visit(p.lists->raw,false);
            if(p.left<0)return;
            i64 contrast=sign*p.contrast;
            if(contrast)visit(contrast>0?p.lists->up:p.lists->down,true);
            read(p.left,sign,true,visit);read(p.right,sign,true,visit);
        }else{
            read(p.left,sign,false,visit);read(p.right,sign,false,visit);
        }
    }
};
