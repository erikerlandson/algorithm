/*******
edit_distance: STL and Boost compatible edit distance functions for C++

Copyright (c) 2013 Erik Erlandson

Author:  Erik Erlandson <erikerlandson@yahoo.com>

Distributed under the Boost Software License, Version 1.0.
See accompanying file LICENSE or copy at
http://www.boost.org/LICENSE_1_0.txt
*******/

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "edit_distance_ut.hpp"

BOOST_AUTO_TEST_SUITE(edit_distance_suite)

#if 0
BOOST_AUTO_TEST_CASE(both_empty) {
    BOOST_CHECK_EQUAL(edit_distance("", ""), 0u);
}

BOOST_AUTO_TEST_CASE(one_empty) {
    BOOST_CHECK_EQUAL(edit_distance("", "a"), 1u);
    BOOST_CHECK_EQUAL(edit_distance("", "ab"), 2u);
    BOOST_CHECK_EQUAL(edit_distance("", "abc"), 3u);

    BOOST_CHECK_EQUAL(edit_distance("a", ""), 1u);
    BOOST_CHECK_EQUAL(edit_distance("ab", ""), 2u);
    BOOST_CHECK_EQUAL(edit_distance("abc", ""), 3u);
}

BOOST_AUTO_TEST_CASE(equal_nonempty) {
    BOOST_CHECK_EQUAL(edit_distance("a", "a"), 0u);
    BOOST_CHECK_EQUAL(edit_distance("ab", "ab"), 0u);
    BOOST_CHECK_EQUAL(edit_distance("abc", "abc"), 0u);
}
#endif

BOOST_AUTO_TEST_CASE(interior_0) {
    BOOST_CHECK_EQUAL(edit_distance("axc", "mxn"), 4u);
    BOOST_CHECK_EQUAL(edit_distance("abxc", "mnxo"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("axbc", "mxno"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("abxc", "mxno"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("axbc", "mnxo"), 6u);
}

BOOST_AUTO_TEST_CASE(interior_1) {
    BOOST_CHECK_EQUAL(edit_distance("a", "b"), 2u);
    BOOST_CHECK_EQUAL(edit_distance("ab", "cd"), 4u);
    BOOST_CHECK_EQUAL(edit_distance("abc", "def"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("abcd", "efgh"), 8u);
    BOOST_CHECK_EQUAL(edit_distance("abcde", "fghij"), 10u);
    BOOST_CHECK_EQUAL(edit_distance("abcdef", "ghijkl"), 12u);
}

BOOST_AUTO_TEST_CASE(interior_2) {
    BOOST_CHECK_EQUAL(edit_distance("a", "bc"), 3u);
    BOOST_CHECK_EQUAL(edit_distance("a", "bcd"), 4u);
    BOOST_CHECK_EQUAL(edit_distance("a", "bcde"), 5u);
    BOOST_CHECK_EQUAL(edit_distance("a", "bcdef"), 6u);

    BOOST_CHECK_EQUAL(edit_distance("ab", "cde"), 5u);
    BOOST_CHECK_EQUAL(edit_distance("ab", "cdef"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("ab", "cdefg"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("ab", "cdefgh"), 8u);

    BOOST_CHECK_EQUAL(edit_distance("abc", "defg"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("abc", "defgh"), 8u);
    BOOST_CHECK_EQUAL(edit_distance("abc", "defghi"), 9u);
    BOOST_CHECK_EQUAL(edit_distance("abc", "defghij"), 10u);
}

BOOST_AUTO_TEST_CASE(interior_2_switch) {
    BOOST_CHECK_EQUAL(edit_distance("ab", "m"), 3u);
    BOOST_CHECK_EQUAL(edit_distance("abc", "m"), 4u);
    BOOST_CHECK_EQUAL(edit_distance("abcd", "m"), 5u);
    BOOST_CHECK_EQUAL(edit_distance("abcde", "m"), 6u);

    BOOST_CHECK_EQUAL(edit_distance("ab", "mn"), 4u);
    BOOST_CHECK_EQUAL(edit_distance("abc", "mn"), 5u);
    BOOST_CHECK_EQUAL(edit_distance("abcd", "mn"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("abcde", "mn"), 7u);

    BOOST_CHECK_EQUAL(edit_distance("ab", "mno"), 5u);
    BOOST_CHECK_EQUAL(edit_distance("abc", "mno"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("abcd", "mno"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("abcde", "mno"), 8u);
}

BOOST_AUTO_TEST_CASE(interior_3) {
    BOOST_CHECK_EQUAL(edit_distance("axbyc", "mxnyo"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("axbyc", "mxnyop"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("axbyc", "mxnyopq"), 8u);
    BOOST_CHECK_EQUAL(edit_distance("axbyc", "mxnyopqr"), 9u);

    BOOST_CHECK_EQUAL(edit_distance("axbyzc", "mxnyzo"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("axbyzc", "mxnyzop"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("axbyzc", "mxnyzopq"), 8u);
    BOOST_CHECK_EQUAL(edit_distance("axbyzc", "mxnyzopqr"), 9u);

    BOOST_CHECK_EQUAL(edit_distance("awxbyzc", "mwxnyzo"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("awxbyzc", "mwxnyzop"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("awxbyzc", "mwxnyzopq"), 8u);
    BOOST_CHECK_EQUAL(edit_distance("awxbyzc", "mwxnyzopqr"), 9u);
}

BOOST_AUTO_TEST_CASE(interior_3_switch) {
    BOOST_CHECK_EQUAL(edit_distance("axbyc", "mxnyo"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("axbycd", "mxnyo"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("axbycde", "mxnyo"), 8u);
    BOOST_CHECK_EQUAL(edit_distance("axbycdef", "mxnyo"), 9u);

    BOOST_CHECK_EQUAL(edit_distance("axbyzc", "mxnyzo"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("axbyzcd", "mxnyzo"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("axbyzcde", "mxnyzo"), 8u);
    BOOST_CHECK_EQUAL(edit_distance("axbyzcdef", "mxnyzo"), 9u);

    BOOST_CHECK_EQUAL(edit_distance("awxbyzc", "mwxnyzo"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("awxbyzcd", "mwxnyzo"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("awxbyzcde", "mwxnyzo"), 8u);
    BOOST_CHECK_EQUAL(edit_distance("awxbyzcdef", "mwxnyzo"), 9u);
}


BOOST_AUTO_TEST_CASE(interior_4) {
    BOOST_CHECK_EQUAL(edit_distance("waaa1bbbx", "yaaabbbzzz"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("waaa12bbbx", "yaaabbbzzz"), 8u);
    BOOST_CHECK_EQUAL(edit_distance("waaa123bbbx", "yaaabbbzzzz"), 10u);
}

BOOST_AUTO_TEST_CASE(interior_4_switch) {
    BOOST_CHECK_EQUAL(edit_distance("waaa1bbbxxx", "yaaabbbz"), 7u);
    BOOST_CHECK_EQUAL(edit_distance("waaa12bbbxxx", "yaaabbbz"), 8u);
}

BOOST_AUTO_TEST_CASE(interior_5) {
    BOOST_CHECK_EQUAL(edit_distance("xxa", "mxx"), 2u);
    BOOST_CHECK_EQUAL(edit_distance("xxab", "mnxx"), 4u);
    BOOST_CHECK_EQUAL(edit_distance("xxabc", "mnoxx"), 6u);
    BOOST_CHECK_EQUAL(edit_distance("xxabcd", "mnopxx"), 8u);
    BOOST_CHECK_EQUAL(edit_distance("xxabc", "mnopxx"), 7u);
}

BOOST_AUTO_TEST_CASE(interior_6) {
    BOOST_CHECK_EQUAL(edit_distance("mxxxxn", "axxxbxx"), 5u);
    BOOST_CHECK_EQUAL(edit_distance("axxxbxx", "mxxxxn"), 5u);
}


#if 0
BOOST_AUTO_TEST_CASE(both_empty_sub) {
    BOOST_CHECK_EQUAL(edit_distance("", "", _substitution=true_type()), (unsigned)0);
}

BOOST_AUTO_TEST_CASE(one_empty_sub) {
    BOOST_CHECK_EQUAL(edit_distance("", "abc", _substitution=true_type()), (unsigned)3);
    BOOST_CHECK_EQUAL(edit_distance("abc", "", _substitution=true_type()), (unsigned)3);
}

BOOST_AUTO_TEST_CASE(length_1_sub) {
    // some boundary conditions for sequence length
    BOOST_CHECK_EQUAL(edit_distance("a", "", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("ab", "a", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("", "a", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("a", "ab", _substitution=true_type()), (unsigned)1);
}

BOOST_AUTO_TEST_CASE(equal_nonempty_sub) {
    BOOST_CHECK_EQUAL(edit_distance("a", "a", _substitution=true_type()), (unsigned)0);
    BOOST_CHECK_EQUAL(edit_distance("ab", "ab", _substitution=true_type()), (unsigned)0);
    BOOST_CHECK_EQUAL(edit_distance("abc", "abc", _substitution=true_type()), (unsigned)0);
}

BOOST_AUTO_TEST_CASE(insertion_sub) {
    // insertion occurs wrt seq2
    BOOST_CHECK_EQUAL(edit_distance("abc", "abcx", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("abc", "abxc", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("abc", "axbc", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("abc", "xabc", _substitution=true_type()), (unsigned)1);

    BOOST_CHECK_EQUAL(edit_distance("abc", "abcxx", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("abc", "abxxc", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("abc", "axxbc", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("abc", "xxabc", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("abc", "axbxc", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("abc", "xabcx", _substitution=true_type()), (unsigned)2);
}

BOOST_AUTO_TEST_CASE(deletion_sub) {
    // deletion occurs wrt seq1
    BOOST_CHECK_EQUAL(edit_distance("abcx", "abc", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("abxc", "abc", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("axbc", "abc", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("xabc", "abc", _substitution=true_type()), (unsigned)1);

    BOOST_CHECK_EQUAL(edit_distance("abcxx", "abc", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("abxxc", "abc", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("axxbc", "abc", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("xxabc", "abc", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("axbxc", "abc", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("xabcx", "abc", _substitution=true_type()), (unsigned)2);
}

BOOST_AUTO_TEST_CASE(substitution_sub) {
    BOOST_CHECK_EQUAL(edit_distance("abc", "axc", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("axc", "abc", _substitution=true_type()), (unsigned)1);
}

BOOST_AUTO_TEST_CASE(sequence_variations) {
    BOOST_CHECK_EQUAL(edit_distance("abc", "axc", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance(ASSTRING("abc"), ASSTRING("axc"), _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance(ASLIST("abc"), ASLIST("axc"), _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance(ASVECTOR("abc"), ASVECTOR("axc"), _substitution=true_type()), (unsigned)1);
}

BOOST_AUTO_TEST_CASE(mixed_sequences) {
    BOOST_CHECK_EQUAL(edit_distance("abc", "bcd", _substitution=true_type()), (unsigned)2);

    BOOST_CHECK_EQUAL(edit_distance("abc", ASSTRING("bcd"), _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("abc", ASLIST("bcd"), _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("abc", ASVECTOR("bcd"), _substitution=true_type()), (unsigned)2);

    BOOST_CHECK_EQUAL(edit_distance(ASSTRING("abc"), "bcd", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance(ASLIST("abc"), "bcd", _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance(ASVECTOR("abc"), "bcd", _substitution=true_type()), (unsigned)2);

    BOOST_CHECK_EQUAL(edit_distance(ASSTRING("abc"), ASLIST("bcd"), _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance(ASVECTOR("abc"), ASLIST("bcd"), _substitution=true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance(ASLIST("abc"), ASVECTOR("bcd"), _substitution=true_type()), (unsigned)2);
}

BOOST_AUTO_TEST_CASE(range_adaptors) {
    BOOST_CHECK_EQUAL(edit_distance("abc", ASLIST("abc") | boost::adaptors::reversed, _substitution=true_type()), (unsigned)2);
}

BOOST_AUTO_TEST_CASE(custom_cost) {
    // make subsitution too expensive to use, so cheapest edit sequence
    // is to delete 'b' and insert 'x'
    BOOST_CHECK_EQUAL(edit_distance("abc", "axc", cost_expensive_sub(), _substitution=true_type()), 2);

    // insertion costs twice as much as deletion: an example of
    // an asymmetric cost function that causes edit distance to be
    // asymmetric
    BOOST_CHECK_EQUAL(edit_distance("aaaa", "aa", cost_expensive_ins(), _substitution=true_type()), 2);
    BOOST_CHECK_EQUAL(edit_distance("aa", "aaaa", cost_expensive_ins(), _substitution=true_type()), 4);
}

BOOST_AUTO_TEST_CASE(undefined_sub) {
    // verify that substitution() method can be undefined when substitution is compile-time disabled:
    BOOST_CHECK_EQUAL(edit_distance("abc", "axc", _cost=undef_sub_cost()), (unsigned)2);
}

BOOST_AUTO_TEST_CASE(custom_equal) {
    BOOST_CHECK_EQUAL(edit_distance("abc", "aBc"), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("abc", "aBc", _equal=case_equal()), (unsigned)0);

    BOOST_CHECK_EQUAL(edit_distance("abc", "aBc", _substitution=true_type()), (unsigned)1);
    BOOST_CHECK_EQUAL(edit_distance("abc", "aBc", _substitution=true_type(), _equal=case_equal()), (unsigned)0);
}

BOOST_AUTO_TEST_CASE(allow_sub_1) {
    // bool arg is run-time check
    BOOST_CHECK_EQUAL(edit_distance("abc", "xyz", _substitution=true), (unsigned)3);
    BOOST_CHECK_EQUAL(edit_distance("abc", "xyz", _substitution=false), (unsigned)6);

    // type arg gives compile-time value: check can be optimized out:
    BOOST_CHECK_EQUAL(edit_distance("abc", "xyz", _substitution=boost::true_type()), (unsigned)3);
    BOOST_CHECK_EQUAL(edit_distance("abc", "xyz", _substitution=boost::false_type()), (unsigned)6);

    BOOST_CHECK_EQUAL(edit_distance("aqc", "xqz", _substitution=boost::true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("aqc", "xqz", _substitution=boost::false_type()), (unsigned)4);

    BOOST_CHECK_EQUAL(edit_distance("aqcr", "xqzr", _substitution=boost::true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("aqcr", "xqzr", _substitution=boost::false_type()), (unsigned)4);

    BOOST_CHECK_EQUAL(edit_distance("raqc", "rxqz", _substitution=boost::true_type()), (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("raqc", "rxqz", _substitution=boost::false_type()), (unsigned)4);
}

#if 0
BOOST_AUTO_TEST_CASE(max_cost_error_1) {
    std::string seq1 = "xx21fxxxxxxxxxxxxxxxxxxxgxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx2b82cfxxxxxxxxxxxxxxxxxxxx";
    std::string seq2 = "x32gc3eaxxxxxxxxxxxxxxxxxedfxxxxxxxxxxxxxxxxxg63bxxxxxxxxxxxxxxxxxxg2xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";

    int d = edit_distance(seq1, seq2);
    BOOST_CHECK_EQUAL(d, 39);

    int dt = edit_distance(seq1, seq2, _max_cost=d-2);
    BOOST_CHECK_GT(dt, d-2); 
}

BOOST_AUTO_TEST_CASE(max_cost_1) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 10000;
    random_localized_deviations(seqdata, N, 100, 5, 10);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            unsigned int d = edit_distance(seqdata[i], seqdata[j]);
            if (d < 2) continue;

            unsigned int dub = seqdata[i].size()+seqdata[j].size();
            unsigned int dt;

            dt = edit_distance(seqdata[i], seqdata[j], _max_cost=d-2);
            BOOST_CHECK_GT(dt, d-2);
            BOOST_CHECK_MESSAGE(dt > d-2, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j] <<"'\n\n");
            BOOST_CHECK_GE(dt, d);
            BOOST_CHECK_LE(dt, dub);

            dt = edit_distance(seqdata[i], seqdata[j], _max_cost=d-1);
            BOOST_CHECK_GT(dt, d-1);
            BOOST_CHECK_GE(dt, d);
            BOOST_CHECK_LE(dt, dub);
            
            dt = edit_distance(seqdata[i], seqdata[j], _max_cost=d);
            BOOST_CHECK_EQUAL(dt, d);

            dt = edit_distance(seqdata[i], seqdata[j], _max_cost=d+1);
            BOOST_CHECK_EQUAL(dt, d);

            dt = edit_distance(seqdata[i], seqdata[j], _max_cost=d+2);
            BOOST_CHECK_EQUAL(dt, d);

            BOOST_CHECK_THROW(edit_distance(seqdata[i], seqdata[j], _max_cost=d-1, _max_cost_exception=true), max_edit_cost_exception);

            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n");
}


BOOST_AUTO_TEST_CASE(max_cost_2) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 200;
    random_localized_deviations(seqdata, N, 100, 5, 10);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            unsigned int d = edit_distance(seqdata[i], seqdata[j], _substitution=true_type());
            if (d < 2) continue;

            unsigned int dub = std::max(seqdata[i].size(), seqdata[j].size());
            unsigned int dt;

            dt = edit_distance(seqdata[i], seqdata[j], _substitution=true_type(), _max_cost=d-2);
            BOOST_CHECK_GT(dt, d-2);
            BOOST_CHECK_GE(dt, d);
            BOOST_CHECK_LE(dt, dub);

            dt = edit_distance(seqdata[i], seqdata[j], _substitution=true_type(), _max_cost=d-1);
            BOOST_CHECK_GT(dt, d-1);
            BOOST_CHECK_GE(dt, d);
            BOOST_CHECK_LE(dt, dub);
            
            dt = edit_distance(seqdata[i], seqdata[j], _substitution=true_type(), _max_cost=d);
            BOOST_CHECK_EQUAL(dt, d);

            dt = edit_distance(seqdata[i], seqdata[j], _substitution=true_type(), _max_cost=d+1);
            BOOST_CHECK_EQUAL(dt, d);

            dt = edit_distance(seqdata[i], seqdata[j], _substitution=true_type(), _max_cost=d+2);
            BOOST_CHECK_EQUAL(dt, d);

            BOOST_CHECK_THROW(edit_distance(seqdata[i], seqdata[j], _substitution=true_type(), _max_cost=d-1, _max_cost_exception=true), max_edit_cost_exception);

            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n");
}

#endif

BOOST_AUTO_TEST_CASE(long_sequences) {
    BOOST_CHECK_EQUAL(edit_distance("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",
                                    "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", _substitution=true_type()),
                      (unsigned)0);
    BOOST_CHECK_EQUAL(edit_distance("xxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",
                                    "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxx", _substitution=true_type()),
                      (unsigned)2);
    BOOST_CHECK_EQUAL(edit_distance("xxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",
                                    "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxx", _substitution=true_type()),
                      (unsigned)4);
    BOOST_CHECK_EQUAL(edit_distance("xxxxxxxxxxxxxxxxxxxxxxxx**xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx**xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",
                                    "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx**xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx**xxxxxxxxxxxxxxxxxxxxxxxxx", _substitution=true_type()),
                      (unsigned)8);
}



BOOST_AUTO_TEST_CASE(myers_sssp_crosscheck_1) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 10;
    random_localized_deviations(seqdata, N, 100000, 5, 50, 100);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            unsigned int d1 = edit_distance(seqdata[i], seqdata[j], _substitution=boost::false_type());
            unsigned int d2 = edit_distance(seqdata[i], seqdata[j], _substitution=boost::false_type(), _cost=unit_cost_test());
            BOOST_CHECK_EQUAL(d1, d2);
            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n");
}


#endif

BOOST_AUTO_TEST_CASE(myers_sssp_crosscheck_2) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 1000;
    random_localized_deviations(seqdata, N, 100, 5, 10);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            unsigned int d1 = edit_distance(seqdata[i], seqdata[j]);
            unsigned int d2 = edit_distance(seqdata[i], seqdata[j], _cost=unit_cost_test());
            BOOST_CHECK_MESSAGE(d1==d2, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j]);
            BOOST_CHECK_EQUAL(d1, d2);
            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n");
}


BOOST_AUTO_TEST_CASE(timing_1) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 100;
    random_localized_deviations(seqdata, N, 1000000, 5, 1000, 100);
    int n = 0;
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            unsigned int d = edit_distance(seqdata[i], seqdata[j]);
            BOOST_CHECK(d <= seqdata[i].size() + seqdata[j].size());
            if (++n >= N) break;
        }
    }
}



#if 0
BOOST_AUTO_TEST_CASE(timing_1_sub) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 100;
    random_localized_deviations(seqdata, N, 100000, 5, 20, 10);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            unsigned int d = edit_distance(seqdata[i], seqdata[j], _substitution=true_type());
            BOOST_CHECK(d <= std::max(seqdata[i].size(),seqdata[j].size()));
            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n" );
}


BOOST_AUTO_TEST_CASE(timing_2_sub) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 100000;
    random_localized_deviations(seqdata, N, 100, 2, 5);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            unsigned int d = edit_distance(seqdata[i], seqdata[j], _substitution=true_type());
            BOOST_CHECK(d <= std::max(seqdata[i].size(),seqdata[j].size()));
            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n" );
}


BOOST_AUTO_TEST_CASE(timing_3_sub) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 1000000;
    random_localized_deviations(seqdata, N, 10, 2, 2);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            unsigned int d = edit_distance(seqdata[i], seqdata[j], _substitution=true_type());
            BOOST_CHECK(d <= std::max(seqdata[i].size(),seqdata[j].size()));
            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n" );
}
#endif


BOOST_AUTO_TEST_SUITE_END()
