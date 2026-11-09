#include "review_base.inc"
#include "drift.hpp"
#include "drift_controls.hpp"
#include "bands.hpp"
#include "drift_split.hpp"
#include <functional>
#include "prefix_class.hpp"

struct QueryWork {DriftWork head;DriftControlWork control;Metrics tail;PrefixWork prefix;};
using Kernel=std::function<i64(const Query&,const vector<Hit>&,QueryWork&)>;

int main(int argc,char** argv) try {
    require(argc==9,"isolated facts bins groups method schedule oracle observations build");
    const std::string method=argv[4];const int bins=std::stoi(argv[2]),groups=std::stoi(argv[3]);
    std::ifstream input(argv[1],std::ios::ate|std::ios::binary);const auto bytes=input.tellg();
    require(bytes>=0 && bytes%sizeof(Fact)==0 && size_t(bytes/sizeof(Fact))<INT_MAX,"invalid input length");
    vector<Fact> facts(size_t(bytes/sizeof(Fact)));input.seekg(0);input.read(reinterpret_cast<char*>(facts.data()),bytes);
    require(bool(input),"truncated facts");
    int last=-1;for(auto f:facts){require(f.time>=last&&f.time<bins&&f.key>=0&&f.key<groups,"invalid fact");last=f.time;}
    i64 retained=0, logical=0;Kernel kernel;
    const auto start=Clock::now();
    if(method.rfind("split",0)==0) {
        const int mode=method=="split-fixed"?2:method=="split-box"?1:method=="split-flat"?0:3;
        require(method=="split"||method=="split-separate"||method=="split-fixed"||method=="split-box"||method=="split-flat","unknown split");
        const bool share=method!="split-separate";
        auto op=std::make_shared<Split>(facts,bins,groups);
        retained=op->bytes();logical=op->bytes(mode);
        kernel=[op,mode,share,method](const Query& q,const vector<Hit>& expected,QueryWork& w) {
            return timed_check([&]{return op->query(q,mode,share,w.head,w.tail);},expected,q.id,method);
        };
    } else if(method=="prefix256" || method=="prefixfull") {
        size_t budget=method=="prefix256" ? size_t(256)<<20 : size_t(8)*bins*groups;
        auto op=std::make_shared<PrefixControl>(facts,bins,groups,budget);
        retained=logical=op->bytes();
        kernel=[op,method](const Query& q,const vector<Hit>& expected,QueryWork& w) {
            return timed_check([&]{return op->query(q,w.prefix);},expected,q.id,method);
        };
    } else if(method=="drift") {
        auto op=std::make_shared<Drift>(facts,bins,groups);retained=logical=op->bytes();
        kernel=[op,method](const Query& q,const vector<Hit>& expected,QueryWork& w) {
            return timed_check([&]{return op->query(q,3,w.head);},expected,q.id,method);
        };
    } else if(method=="band-raw") {
        auto op=std::make_shared<Bands>(facts,bins,groups,"band-raw");retained=logical=op->bytes();
        kernel=[op,method](const Query& q,const vector<Hit>& expected,QueryWork& w) {
            return timed_check([&]{return op->query(q,w.tail);},expected,q.id,method);
        };
    } else {
        require(method=="ctrl1"||method=="ctrl2","unknown control");const int mode=method.back()-'0';
        auto op=std::make_shared<DriftControls>(facts,bins,groups,mode);retained=logical=op->bytes();
        kernel=[op,mode,method](const Query& q,const vector<Hit>& expected,QueryWork& w) {
            return timed_check([&]{return op->query(q,mode,w.control);},expected,q.id,method);
        };
    }
    const auto build_ns=std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-start).count();
    vector<Fact>().swap(facts);
    PROCESS_MEMORY_COUNTERS memory{};memory.cb=sizeof(memory);
    require(GetProcessMemoryInfo(GetCurrentProcess(),&memory,sizeof(memory))!=0,"cannot measure process memory");
    std::ofstream build(argv[8]);build<<"method,build_ns,retained_bytes,mode_bytes,postbuild_rss_bytes,construction_peak_rss_bytes\n"
        <<method<<','<<build_ns<<','<<retained<<','<<logical<<','<<memory.WorkingSetSize<<','<<memory.PeakWorkingSetSize<<'\n';
    auto oracle=read_oracle(argv[6]);const vector<Hit> empty;
    std::ifstream schedule(argv[5]);std::ofstream obs(argv[7]);std::string line;std::getline(schedule,line);
    obs<<"round,step,qid,method,ns,head_blocks,box_skips,joint_skips,head_scores,cold_scores,cold_cells,tail_verified,tail_ids,tail_passes,control_block,control_key,control_scores,sorted_reads,verified,residual_facts,full\n";
    while(std::getline(schedule,line)) {
        std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);int round,step;Query q;
        require(bool(s>>round>>step>>q.id>>q.a>>q.b>>q.c>>q.d>>q.k),"bad query record");
        require(q.a>=0&&q.a<=q.b&&q.b<=bins&&q.c>=0&&q.c<=q.d&&q.d<=bins&&q.k>0,"bad bounds");
        auto it=oracle.find(q.id);const auto& expected=it==oracle.end()?empty:it->second;
        QueryWork w{};const auto ns=kernel(q,expected,w);
        obs<<round<<','<<step<<','<<q.id<<','<<method<<','<<ns<<','<<w.head.seen<<','<<w.head.box<<','<<w.head.pair<<','<<w.head.scored<<','<<w.head.cold<<','<<w.head.cells<<','<<w.tail.verified<<','<<w.tail.ids<<','<<w.tail.passes<<','<<w.control.block_skips<<','<<w.control.key_skips<<','<<w.control.hot_scores+w.control.cold_scores<<','<<w.prefix.reads<<','<<w.prefix.verified<<','<<w.prefix.residual_facts<<','<<w.prefix.fallback<<'\n';
    }
    build.close();obs.close();require(bool(build)&&bool(obs),"writing failed");
} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 2; }
