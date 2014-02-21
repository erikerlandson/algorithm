/*******
edit_distance: STL and Boost compatible edit distance functions for C++

Copyright (c) 2013 Erik Erlandson

Author:  Erik Erlandson <erikerlandson@yahoo.com>

Distributed under the Boost Software License, Version 1.0.
See accompanying file LICENSE or copy at
http://www.boost.org/LICENSE_1_0.txt
*******/


#if !defined(BOOST_ALGORITHM_SEQUENCE_DETAIL_EDIT_DISTANCE_HPP)
#define BOOST_ALGORITHM_SEQUENCE_DETAIL_EDIT_DISTANCE_HPP

// this requires me to link against boost system lib, which disappoints me,
// since it prevents my algorithms from being pure-template.
// maybe figure out how to use a pure-template pool allocator later
#include <boost/pool/object_pool.hpp>

#include <boost/heap/skew_heap.hpp>

#include <boost/range/metafunctions.hpp>

#include <boost/algorithm/sequence/detail/edit_distance_types.hpp>
#include <boost/algorithm/sequence/detail/edit_distance_script.hpp>

namespace boost {
namespace algorithm {
namespace sequence {
namespace detail {

using boost::distance;
using boost::range_iterator;

using boost::make_tuple;
using boost::enable_if;
using boost::is_same;
using boost::mpl::and_;
using std::iterator_traits;
using std::random_access_iterator_tag;

template <typename ForwardRange1, typename ForwardRange2, typename Output, typename Cost, typename Equal, typename AllowSub, typename MaxCost>
struct edit_cost_struct<ForwardRange1, ForwardRange2, Output, Cost, Equal, AllowSub, MaxCost,
                        typename enable_if<and_<is_same<Output, none>,
                                                not_<and_<range_category<ForwardRange1, ForwardRange2, random_access_iterator_tag>,
                                                          is_same<Cost, unit_cost>,
                                                          is_same<AllowSub, boost::false_type> > > > >::type> {
typedef typename cost_type<Cost, typename boost::range_value<ForwardRange1>::type>::type cost_t;
typedef typename range_iterator<ForwardRange1 const>::type itr1_t;
typedef typename range_iterator<ForwardRange2 const>::type itr2_t;
typedef path_head<itr1_t, itr2_t, cost_t> head_t;
typedef typename head_t::pos1_type pos1_t;
typedef typename head_t::pos2_type pos2_t;

cost_t max_cost_fallback(max_cost_checker<MaxCost, cost_t, head_t>& max_cost_check, bool max_cost_exception, const itr1_t end1, const itr2_t end2, const Cost& cost, const Equal& equal, sub_checker<AllowSub, Cost, cost_t, int> const& allow_sub) const {
    if (max_cost_exception) throw max_edit_cost_exception();

    head_t* h;
    max_cost_check.get(h);

    pos1_t j1 = h->pos1;
    pos2_t j2 = h->pos2;
    cost_t C = h->cost;
    while (true) {
        if (j1 == end1) {
            if (j2 == end2) {
                return C;
            } else {
                C += cost.insertion(*j2);
                ++j2;
            }
        } else {
            if (j2 == end2) {
                C += cost.deletion(*j1);
                ++j1;
            } else {
                C += (equal(*j1, *j2)) ? 0 
                                       : ((allow_sub()) ? std::min(allow_sub.substitution(cost, *j1, *j2), (cost.deletion(*j1)+cost.insertion(*j2))) 
                                                        : (cost.deletion(*j1)+cost.insertion(*j2))) ;
                ++j1;  ++j2;
            }
        }
    }
    return C;
}

// Default is generic edit distance algorithm based on a Dijkstra Single Source Shortest Path approach
typename cost_type<Cost, typename boost::range_value<ForwardRange1>::type>::type
operator()(ForwardRange1 const& seq1, ForwardRange2 const& seq2, none&, const Cost& cost, const Equal& equal, const AllowSub& allowsub, const MaxCost& max_cost, const bool max_cost_exception) const {

    head_t* const hnull = static_cast<head_t*>(NULL);

    const itr1_t end1 = boost::end(seq1);
    const itr2_t end2 = boost::end(seq2);
    pos1_t beg1;  beg1.beg(boost::begin(seq1));
    pos2_t beg2;  beg2.beg(boost::begin(seq2));

    // pool allocator for path nodes
    boost::object_pool<head_t> pool;

    // priority queue for path nodes
    boost::heap::skew_heap<head_t*, boost::heap::compare<heap_lessthan<pos1_t, pos2_t> > > heap(heap_lessthan<pos1_t, pos2_t>(beg1, beg2));

    sub_checker<AllowSub, Cost, cost_t, int> allow_sub(allowsub);

    max_cost_checker<MaxCost, cost_t, head_t> max_cost_check(max_cost, beg1, beg2);

    // keep track of nodes in the edit graph that have been visited
    typedef boost::unordered_set<head_t*, visited_hash<pos1_t,pos2_t>, visited_equal> visited_t;
    visited_t visited(31, visited_hash<pos1_t,pos2_t>(beg1,beg2));

    // kick off graph path frontier with initial node:
    heap.push(construct(pool, visited, beg1, beg2, cost_t(0)));

    // update frontier from least-cost node at each iteration, until we hit sequence end
    while (true) {
        head_t* h = heap.top();
        heap.pop();

        if (max_cost_check(h->cost)) {
            return max_cost_fallback(max_cost_check, max_cost_exception, end1, end2, cost, equal, allow_sub);
        }
        max_cost_check.update(h);

        if (h->pos1 == end1) {
            // if we are at end of both sequences, then we have our final cost: 
            if (h->pos2 == end2) return h->cost;
            // sequence 1 is at end, so only consider insertion from seq2
            pos2_t p2 = h->pos2;
            head_t* t = construct(pool, visited, h->pos1, ++p2, h->cost + cost.insertion(*(h->pos2)));
            if (t != hnull) heap.push(t);
        } else if (h->pos2 == end2) {
            // sequence 2 is at end, so only consider deletion from seq1
            pos1_t p1 = h->pos1;
            head_t* t = construct(pool, visited, ++p1, h->pos2, h->cost + cost.deletion(*(h->pos1)));
            if (t != hnull) heap.push(t);
        } else {
            // interior of both sequences: consider insertion deletion and sub/eql:
            pos1_t p1 = h->pos1;  ++p1;
            pos1_t p1p = h->pos1;
            pos2_t p2 = h->pos2;  ++p2;
            pos2_t p2p = h->pos2;
            while (true) {
                const bool eq = equal(*p1p, *p2p);
                if (!eq  ||  p1 == end1  ||  p2 == end2) {
                    head_t* t;
                    if (allow_sub() || eq) {
                        t = construct(pool, visited, p1, p2, h->cost + ((eq) ? 0 : allow_sub.substitution(cost, *p1p, *p2p)));
                        if (t != hnull) heap.push(t);
                    }
                    t = construct(pool, visited, p1p, p2, h->cost + cost.insertion(*p2p));
                    if (t != hnull) heap.push(t);
                    t = construct(pool, visited, p1, p2p, h->cost + cost.deletion(*p1p));
                    if (t != hnull) heap.push(t);
                    break;
                }
                ++p1;  ++p2;  ++p1p;  ++p2p;
            }
        }
    }

    // control should not reach here
    BOOST_ASSERT(false);
    return 0;
}

}; // edit_cost_struct


template <typename Range1, typename Range2, typename Output, typename Equal, typename MaxCost>
struct edit_cost_struct<Range1, Range2, Output, unit_cost, Equal, boost::false_type, MaxCost,
                        typename enable_if<and_<is_same<Output, none>,
                                                range_category<Range1, Range2, random_access_iterator_tag> > >::type> {

typedef typename range_iterator<Range1 const>::type itr1_t;
typedef typename range_iterator<Range2 const>::type itr2_t;

typedef std::vector<int>::difference_type diff_type;

typedef typename std::vector<diff_type>::iterator itrv_t;
typedef max_cost_checker_myers<MaxCost, diff_type, diff_type> max_cost_type;

std::string dump(const itr1_t& S1, const diff_type& len1) const {
    std::string r;
    for (int j = 0; j < len1;  ++j) r += S1[j];
    return r;
}

void dump(const string& lab, diff_type k, const itr1_t& S1, const diff_type& len1, const itr2_t& S2, const diff_type& len2, diff_type j1, diff_type j2, diff_type t1, diff_type t2) const {
    if (len2 >= len1) {
        std::cout << "    " << make_tuple(lab, k, j1, j2, (j1<len1)?S1[j1]:'_', (j2<len2)?S2[j2]:'_', t1, t2, (t1 < len1)?S1[t1]:'_', (t2 < len2)?S2[t2]:'_') << std::endl;
    } else {
        std::cout << "    " << make_tuple(lab, k, j1, j2, (j2<len1)?S1[j2]:'_', (j1<len2)?S2[j1]:'_', t1, t2, (t2 < len1)?S1[t2]:'_', (t1 < len2)?S2[t1]:'_') << std::endl;
    }
}

template <typename Vec, typename Itr> 
inline void expand(Vec& V_data, Itr& Vf, Itr& Vr, diff_type& R, const diff_type& P, const diff_type& delta, const diff_type& L) const {
    diff_type Rp = R + (R>>1);
    V_data.resize(2*(1 + delta + 2*Rp));
    Vf = V_data.begin() + R;
    Vr = V_data.begin() + (1 + delta + 2*R) + R;

    Itr Vp = V_data.begin() + (1 + delta + 2*Rp) + Rp;
    for (diff_type j=P+delta;  j >= -P;  --j) Vp[j] = Vr[j];
    Vr = Vp;

    Vp = V_data.begin() + Rp;
    for (diff_type j=P+delta;  j >= -P;  --j) Vp[j] = Vf[j];
    Vf = Vp;

    R = Rp;

    for (diff_type j=1;  j<=(R-P);  ++j) {
        Vf[-P-j] = Vf[P+delta+j] = -1;
        Vr[-P-j] = Vr[P+delta+j] = L;
    }
}


diff_type max_cost_fallback(max_cost_type& max_cost_check, const bool max_cost_exception, const Equal& equal,
                            const itr1_t& S1, const diff_type& L1, const itr2_t& S2, const diff_type& L2,
                            const itrv_t& Vf, const itrv_t& Vr, const diff_type& delta, const diff_type& D) const {
    if (max_cost_exception) throw max_edit_cost_exception();

    for (diff_type k = -D;  k <= D;  k += 2) {
        max_cost_check.update(k, Vf, Vr, delta, L1, L2, D);
    }

    diff_type r1b=0, r2b=0, r1e=0, r2e=0;

    diff_type C = 0;
    diff_type k = 0;
    remainder::kind kind;
    max_cost_check.get(k, kind);
    switch (kind) {
    case remainder::forward: {
            r1b = Vf[k];
            r2b = r1b-k;
            r1e = L1;
            r2e = L2;
            C = D;
        }; break;

        case remainder::reverse: {
            r1b = 0;
            r2b = 0;
            r1e = Vr[k];
            r2e = r1e-k;
            C = D;
        }; break;

        case remainder::bidirectional: {
            r1b = Vf[k];
            r2b = r1b-k;
            r1e = Vr[k];
            r2e = r1e-k;
            C = 2*D;
        }; break;

        default: BOOST_ASSERT(false);
    }

    // this is the part we bailed on due to hitting the maximum
    diff_type j1 = r1b;
    diff_type j2 = r2b;
    while (true) {
        if (j1 >= r1e) {
            if (j2 >= r2e) {
                break;
            } else {
                C += 1;
                ++j2;
            }
        } else {
            if (j2 >= r2e) {
                C += 1;
                ++j1;
            } else {
                if (!equal(S1[j1], S2[j2])) C += 2;
                ++j1;
                ++j2;
            }
        }
    }

    return C;
}

// If we are using unit cost for ins/del, with no substitution,
// and if our sequences support random-access,
// *then* we can invoke the efficient and elegant Myers algorithm:
//     An O(ND) Difference Algorithm and its Variations
//     by Eugene W. Myers
//     Dept of Computer Science, University of Arizona
typename cost_type<unit_cost, typename boost::range_value<Range1>::type>::type
operator()(Range1 const& seq1_, Range2 const& seq2_, none&, const unit_cost&, const Equal& equal, const boost::false_type&, const MaxCost& max_cost, const bool max_cost_exception) const {
    itr1_t seq1 = boost::begin(seq1_);
    itr2_t seq2 = boost::begin(seq2_);
    diff_type len1 = distance(seq1_);
    diff_type len2 = distance(seq2_);

//    std::cout << std::endl;
//    std::cout << "seq1= " << dump(seq1, len1) << std::endl;
//    std::cout << "seq2= " << dump(seq2, len2) << std::endl;

    // identify any equal suffix and/or prefix
    diff_type eqb = 0;
    for (;  eqb < std::min(len1, len2);  ++eqb) if (!equal(seq1[eqb],seq2[eqb])) break;
    diff_type eqe = len1-1;
    for (diff_type j2 = len2-1;  eqe > eqb && j2 > eqb;  --eqe,--j2) if (!equal(seq1[eqe],seq2[j2])) break;
    eqe = len1-1-eqe;

    // sub-strings with equal suffix and/or prefix stripped
    const itr1_t S1 = seq1 + eqb;
    const diff_type L1 = len1-(eqb+eqe);
    const itr2_t S2 = seq2 + eqb;
    const diff_type L2 = len2-(eqb+eqe);

    // either or both strings are empty:
    if (L1 <= 0) return L2;
    if (L2 <= 0) return L1;

//    std::cout << "S1= " << dump(S1, L1) << std::endl;
//    std::cout << "S2= " << dump(S2, L2) << std::endl;

    const diff_type delta = (L2 >= L1) ? (L2-L1) : (L1-L2);

    diff_type R = 5;
    std::vector<diff_type> V_data(2*(1 + delta + 2*R));
    itrv_t Vf = V_data.begin() + R;
    itrv_t Vr = V_data.begin() + (1 + delta + 2*R) + R;
    for (diff_type k = -R;  k <= delta+R;  ++k) {
        Vf[k] = -1;
        Vr[k] = 1+std::max(L1, L2);
    }

    max_cost_type max_cost_check(max_cost);

    // initialize this with the maximum possible distance:
    diff_type Dbest = L1+L2;

    diff_type P = 0;
    while (true) {
        // the minimum possible distance for the current P value
        diff_type Dmin = 2*(P-1) + delta;

        // if the minimum possible distance is >= our best-known distance, we can halt
        if (Dmin >= Dbest) return Dbest;

        diff_type bound = std::min(delta, ((Dbest-delta)/2)-(2*(P-1)));

        // advance forward diagonals
        for (diff_type ku = -P, kd = P+delta;  ku <= bound;  ++ku) {
            diff_type j2 = std::max(1+Vf[ku-1], Vf[ku+1]);
            diff_type j1 = j2-ku;

            if (j2 >= Vr[ku]) {
                diff_type vf = (ku>delta) ? (P + delta - ku) : P;
                diff_type vr = (ku<0) ? (P-1 + ku) : P-1;
                Dbest = std::min(Dbest, 2*(vf+vr)+delta);
                break;
            }

            if (L2 >= L1) {
                while (j1 < L1  &&  j2 < L2  &&  equal(S1[j1], S2[j2])) { ++j1;  ++j2; }
            } else {
                while (j1 < L2  &&  j2 < L1  &&  equal(S1[j2], S2[j1])) { ++j1;  ++j2; }
            }

            Vf[ku] = j2;

            if (kd <= delta) continue;

            j2 = std::max(1+Vf[kd-1], Vf[kd+1]);
            j1 = j2-kd;

            if (j2 >= Vr[kd]) {
                diff_type vf = (kd>delta) ? (P + delta - kd) : P;
                diff_type vr = (kd<0) ? (P-1 + kd) : P-1;
                Dbest = std::min(Dbest, 2*(vf+vr)+delta);
                break;
            }

            if (L2 >= L1) {
                while (j1 < L1  &&  j2 < L2  &&  equal(S1[j1], S2[j2])) { ++j1;  ++j2; }
            } else {
                while (j1 < L2  &&  j2 < L1  &&  equal(S1[j2], S2[j1])) { ++j1;  ++j2; }
            }

            Vf[kd] = j2;
            --kd;
        }

        bound = std::max(diff_type(0), ((delta-Dbest)/2)+delta+(2*P));

        // advance reverse-path diagonals:
        for (diff_type kd=P+delta, ku=-P;  kd >= bound;  --kd) {
            diff_type j2 = std::min(Vr[kd-1], Vr[kd+1]-1);
            diff_type j1 = j2-kd;

            if (j2 <= Vf[kd]) {
                diff_type vf = (kd>delta) ? (P + delta - kd) : P;
                diff_type vr = (kd<0) ? (P + kd) : P;
                Dbest = std::min(Dbest, 2*(vf+vr)+delta);
                break;
            }

            if (L2 >= L1) {
                while (j1 > 0  &&  j2 > 0  &&  equal(S1[j1-1], S2[j2-1])) { --j1;  --j2; }
            } else {
                while (j1 > 0  &&  j2 > 0  &&  equal(S1[j2-1], S2[j1-1])) { --j1;  --j2; }
            }

            Vr[kd] = j2;

            if (ku >= 0) continue;

            j2 = std::min(Vr[ku-1], Vr[ku+1]-1);
            j1 = j2-ku;

            if (j2 <= Vf[ku]) {
                diff_type vf = (ku>delta) ? (P + delta - ku) : P;
                diff_type vr = (ku<0) ? (P + ku) : P;
                Dbest = std::min(Dbest, 2*(vf+vr)+delta);
                break;
            }

            if (L2 >= L1) {
                while (j1 > 0  &&  j2 > 0  &&  equal(S1[j1-1], S2[j2-1])) { --j1;  --j2; }
            } else {
                while (j1 > 0  &&  j2 > 0  &&  equal(S1[j2-1], S2[j1-1])) { --j1;  --j2; }
            }

            Vr[ku] = j2;
            ++ku;
        }

#if 0
        if (max_cost_check((delta_even) ? (2*D+2) : (2*D+1))) {
            return max_cost_fallback(max_cost_check, max_cost_exception, equal,
                                     S1, L1, S2, L2,
                                     Vf, Vr, delta, D);
        }
#endif

        // expand the working vector as needed
        if (1+P >= R) expand(V_data, Vf, Vr, R, P, delta, 1+std::max(L1, L2));
        ++P;
    }

    // control should not reach here
    BOOST_ASSERT(false);
    return 0;
}

}; // edit_cost_struct


template <typename Range1, typename Range2, typename Output, typename Cost, typename Equal, typename AllowSub, typename MaxCost>
inline
typename cost_type<Cost, typename boost::range_value<Range1>::type>::type
edit_cost_impl(Range1 const& seq1, Range2 const& seq2, Output& output, const Cost& cost, const Equal& equal, const AllowSub& allow_sub, const MaxCost& max_cost, const bool max_cost_exception) {
    // specialize the most appropriate implementation for the given parameters
    return edit_cost_struct<Range1, Range2, Output, Cost, Equal, AllowSub, MaxCost>()(seq1, seq2, output, cost, equal, allow_sub, max_cost, max_cost_exception);
}


}}}}

#endif
