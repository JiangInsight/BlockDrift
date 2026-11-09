#include "cover.hpp"
#include <iomanip>
#include <sstream>
#include <windows.h>
#include <psapi.h>

struct PrefixColumn {
    vector<int> count, order;
};
struct PrefixWork {
    i64 reads = 0, residual_facts = 0, verified = 0;
    bool fallback = false;
};
struct Cursor {
    const PrefixColumn* list;
    int weight, lo, hi;
};

class PrefixControl {
    int bins, groups, epoch = 0;
    vector<int> keys, offsets, seen, owner;
    vector<i64> adjustment;
    vector<PrefixColumn> columns;

    void snapshot(const vector<int>& cumulative, vector<int>& scratch) {
        PrefixColumn col{cumulative, vector<int>(groups)};
        int maximum = *std::max_element(cumulative.begin(), cumulative.end());


        if (i64(maximum) <= std::max<i64>(65536, 4LL * groups)) {
            scratch.assign(size_t(maximum) + 1, 0);
            for (int value : cumulative) ++scratch[value];
            int next = 0;
            for (int& count : scratch) { int prior = count; count = next; next += prior; }
            for (int k = 0; k < groups; ++k) col.order[scratch[cumulative[k]]++] = k;
        } else {
            std::iota(col.order.begin(), col.order.end(), 0);
            std::sort(col.order.begin(), col.order.end(), [&](int a, int b) {
                return cumulative[a] != cumulative[b] ? cumulative[a] < cumulative[b] : a < b;
            });
        }
        columns.push_back(std::move(col));
    }
public:
    int stride;
    PrefixControl(const vector<Fact>& facts, int t, int g, size_t budget):
        bins(t), groups(g), keys(facts.size()), offsets(t+1), seen(g), owner(t+1,-1), adjustment(g) {
        size_t room = budget / (sizeof(int) * 2ULL * g);
        require(room >= 1, "snapshot budget cannot hold a single full column");
        int max_cols = int(std::min<size_t>(room, t));
        stride = (t + max_cols - 1) / max_cols;
        columns.reserve((t + stride - 1) / stride);
        vector<int> cumulative(g), scratch;
        size_t i = 0;
        for (int end = 1; end <= t; ++end) {
            while (i < facts.size() && facts[i].time == end-1) {
                int k = facts[i].key;
                ++cumulative[k]; keys[i++] = k;
            }
            offsets[end] = int(i);
            owner[end] = owner[end-1];
            if (end % stride == 0 || end == t) {
                owner[end] = int(columns.size());
                snapshot(cumulative, scratch);
            }
        }
    }
    size_t bytes() const {
        size_t integers = keys.capacity()+offsets.capacity()+seen.capacity()+owner.capacity();
        for (const auto& c : columns) integers += c.count.capacity()+c.order.capacity();
        return sizeof(*this)+integers*sizeof(int)+adjustment.capacity()*sizeof(i64)+columns.capacity()*sizeof(PrefixColumn);
    }
    size_t snapshots() const { return columns.size(); }
    vector<Hit> query(const Query& q, PrefixWork& work) {
        ++epoch;
        std::array<std::pair<int,int>,4> endpoints{{{q.a,1},{q.b,-1},{q.c,-1},{q.d,1}}};
        std::sort(endpoints.begin(),endpoints.end());
        vector<std::pair<int,int>> weights;
        vector<std::pair<int,int>> cuts;
        for (size_t a = 0; a < endpoints.size();) {
            int t = endpoints[a].first, sign = 0;
            while (a < endpoints.size() && endpoints[a].first == t) sign += endpoints[a++].second;
            if (!sign) continue;
            int col = owner[t];
            int anchor = t == bins ? bins : (t / stride) * stride;
            if (col != -1) {
                if (!weights.empty() && weights.back().first == col) weights.back().second += sign;
                else weights.emplace_back(col, sign);
            }
            if (anchor != t) { cuts.emplace_back(anchor, sign); cuts.emplace_back(t,-sign); }
        }

        std::sort(cuts.begin(), cuts.end());
        vector<int> residual;
        int weight = 0, prev = 0;
        for (auto [end, change] : cuts) {
            if (weight) for (int i = offsets[prev]; i < offsets[end]; ++i) {
                int k = keys[i];
                if (seen[k] != epoch) { seen[k] = epoch; residual.push_back(k); }
                adjustment[k] += weight; ++work.residual_facts;
            }
            weight += change; prev = end;
        }
        vector<Cursor> cursors;
        for (auto [col, coefficient] : weights) if (coefficient)
            cursors.push_back({&columns[col],coefficient,0,groups-1});
        Top top(q.k);
        auto score = [&](int k) {
            i64 d = adjustment[k];
            for (const Cursor& c : cursors) d += i64(c.weight)*c.list->count[k];
            return int(d);
        };
        auto verify = [&](int k) {
            if (seen[k] == epoch) return;
            seen[k] = epoch; ++work.verified; top.add(k, score(k));
        };
        for (int k : residual) {
            seen[k] = 0;
            if (adjustment[k] != 0) verify(k);
        }
        while (!cursors.empty()) {
            i64 upper = 0, lower = 0;
            bool exhausted = false;
            for (const Cursor& c : cursors) {
                if (c.lo > c.hi) { exhausted = true; break; }
                int lo = c.list->count[c.list->order[c.lo]], hi = c.list->count[c.list->order[c.hi]];
                i64 x = i64(c.weight)*lo, y = i64(c.weight)*hi;
                lower += std::min(x,y); upper += std::max(x,y);
            }
            i64 bound = std::max(upper,-lower);
            if (exhausted || !bound || bound < top.threshold()) break;

            if (work.reads >= groups) {
                work.fallback = true;
                for (int k=0;k<groups;++k) verify(k);
                break;
            }
            for (Cursor& c : cursors) {
                verify(c.list->order[c.lo++]); ++work.reads;
                if (c.lo <= c.hi) { verify(c.list->order[c.hi--]); ++work.reads; }
            }
        }
        vector<Hit> hits=top.finish(score);
        for (int k : residual) adjustment[k] = 0;
        return hits;
    }
};
