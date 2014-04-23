/*******
edit_distance: STL and Boost compatible edit distance functions for C++

Copyright (c) 2013 Erik Erlandson

Author:  Erik Erlandson <erikerlandson@yahoo.com>

Distributed under the Boost Software License, Version 1.0.
See accompanying file LICENSE or copy at
http://www.boost.org/LICENSE_1_0.txt
*******/

#include <boost/test/unit_test.hpp>

#include "edit_distance_ut.hpp"

BOOST_AUTO_TEST_SUITE(edit_alignment_suite)

BOOST_AUTO_TEST_CASE(both_empty) {
    CHECK_EDIT_ALIGNMENT_ARG("", "", _substitution=true_type(), 0);
}

BOOST_AUTO_TEST_CASE(insertion) {
    CHECK_EDIT_ALIGNMENT_ARG("", "a", _substitution=true_type(), 1);
    CHECK_EDIT_ALIGNMENT_ARG("", "aa", _substitution=true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG("", "aaa", _substitution=true_type(), 3);
}

BOOST_AUTO_TEST_CASE(deletion) {
    CHECK_EDIT_ALIGNMENT_ARG("a", "", _substitution=true_type(), 1);
    CHECK_EDIT_ALIGNMENT_ARG("aa", "", _substitution=true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG("aaa", "", _substitution=true_type(), 3);
}

BOOST_AUTO_TEST_CASE(substitution) {
    CHECK_EDIT_ALIGNMENT_ARG("a", "x", _substitution=true_type(), 1);
    CHECK_EDIT_ALIGNMENT_ARG("ab", "xy", _substitution=true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG("abc", "xyz", _substitution=true_type(), 3);
}

BOOST_AUTO_TEST_CASE(substitution_equal) {
    CHECK_EDIT_ALIGNMENT_ARG("a", "a", _substitution=true_type(), 0);
    CHECK_EDIT_ALIGNMENT_ARG("aa", "aa", _substitution=true_type(), 0);
    CHECK_EDIT_ALIGNMENT_ARG("aaa", "aaa", _substitution=true_type(), 0);
}

BOOST_AUTO_TEST_CASE(sequence_variations) {
    CHECK_EDIT_ALIGNMENT_ARG("abc", "axc", _substitution=true_type(), 1);
    CHECK_EDIT_ALIGNMENT_ARG(ASSTRING("abc"), ASSTRING("axc"), _substitution=true_type(), 1);
    CHECK_EDIT_ALIGNMENT_ARG(ASLIST("abc"), ASLIST("axc"), _substitution=true_type(), 1);
    CHECK_EDIT_ALIGNMENT_ARG(ASVECTOR("abc"), ASVECTOR("axc"), _substitution=true_type(), 1);
}

BOOST_AUTO_TEST_CASE(mixed_sequences) {
    CHECK_EDIT_ALIGNMENT_ARG("abc", "bcd", _substitution=true_type(), 2);

    CHECK_EDIT_ALIGNMENT_ARG("abc", ASSTRING("bcd"), _substitution=true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG("abc", ASLIST("bcd"), _substitution=true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG("abc", ASVECTOR("bcd"), _substitution=true_type(), 2);

    CHECK_EDIT_ALIGNMENT_ARG(ASSTRING("abc"), "bcd", _substitution=true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG(ASLIST("abc"), "bcd", _substitution=true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG(ASVECTOR("abc"), "bcd", _substitution=true_type(), 2);

    CHECK_EDIT_ALIGNMENT_ARG(ASSTRING("abc"), ASLIST("bcd"), _substitution=true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG(ASVECTOR("abc"), ASLIST("bcd"), _substitution=true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG(ASLIST("abc"), ASVECTOR("bcd"), _substitution=true_type(), 2);
}

BOOST_AUTO_TEST_CASE(range_adaptors) {
    CHECK_EDIT_ALIGNMENT_ARG("abc", ASLIST("abc") | boost::adaptors::reversed, _substitution=true_type(), 2);
}

BOOST_AUTO_TEST_CASE(mixed_ops) {
    CHECK_EDIT_ALIGNMENT_ARG("abcd", "bCde", _substitution=true_type(), 3);
}

BOOST_AUTO_TEST_CASE(custom_cost) {
    CHECK_EDIT_ALIGNMENT_2ARG("abcd", "bCde", _substitution=true_type(), _cost = cost_expensive_sub(), 4);
    CHECK_EDIT_ALIGNMENT_2ARG("abcd", "aBCd", _substitution=true_type(), _cost = cost_expensive_sub(), 4);

    CHECK_EDIT_ALIGNMENT_2ARG("aa", "axax", _substitution=true_type(), _cost = cost_expensive_ins(), 4);
    CHECK_EDIT_ALIGNMENT_2ARG("axax", "aa", _substitution=true_type(), _cost = cost_expensive_ins(), 2);

    CHECK_EDIT_ALIGNMENT_2ARG("aa", "axax", _substitution=true_type(), _cost = cost_expensive_ins(), 4);
    CHECK_EDIT_ALIGNMENT_2ARG("axax", "aa", _substitution=true_type(), _cost = cost_expensive_ins(), 2);
}

BOOST_AUTO_TEST_CASE(undefined_sub) {
    // verify that substitution() and substitution() can be undefined when substitution is compile-time disabled:
    const std::string seq1 = "abc";
    const std::string seq2 = "axc";

    {
        // test sssp implementation
        undef_sub_output<char, unsigned> ob(ASVECTOR(seq1), ASVECTOR(seq2));
        BOOST_CHECK_EQUAL(edit_distance(seq1, seq2, _script=ob, _cost=undef_sub_cost()), unsigned(2));
        ob.finalize(2);
        BOOST_CHECK(ob.correct);
    }

    {
        // test invoking myers specialization
        undef_sub_output<char, unsigned> ob(ASVECTOR(seq1), ASVECTOR(seq2));
        BOOST_CHECK_EQUAL(edit_distance(seq1, seq2, _script=ob), unsigned(2));
        ob.finalize(2);
        BOOST_CHECK(ob.correct);
    }
}


BOOST_AUTO_TEST_CASE(custom_equal) {
    {
        const std::string seq1 = "abc";
        const std::string seq2 = "aBc";
        const unsigned dist = 2;
        output_check_script<char, unsigned> ob(ASVECTOR(seq1), ASVECTOR(seq2));
        unsigned d = edit_distance(seq1, seq2, _script=ob);
        ob.finalize(dist);
        BOOST_CHECK_MESSAGE(ob.correct, "incorrect edit script: '" << ob.ss.str() << "'  seq1='" << seq1 << "'  seq2='" << seq2 << "'");
        BOOST_CHECK_MESSAGE(d == dist, "incorrect edit distance " << d << "(expected " << dist << ")  seq1='" << seq1 << "' seq2='" << seq2 << "'  script='" << ob.ss.str() <<"'");
    }
    {
        const std::string seq1 = "abc";
        const std::string seq2 = "aBc";
        const unsigned dist = 0;
        output_check_script<char, unsigned, case_equal> ob(ASVECTOR(seq1), ASVECTOR(seq2));
        unsigned d = edit_distance(seq1, seq2, _script=ob, _equal=case_equal());
        ob.finalize(dist);
        BOOST_CHECK_MESSAGE(ob.correct, "incorrect edit script: '" << ob.ss.str() << "'  seq1='" << seq1 << "'  seq2='" << seq2 << "'");
        BOOST_CHECK_MESSAGE(d == dist, "incorrect edit distance " << d << "(expected " << dist << ")  seq1='" << seq1 << "' seq2='" << seq2 << "'  script='" << ob.ss.str() <<"'");
    }
    {
        const std::string seq1 = "abc";
        const std::string seq2 = "aBc";
        const unsigned dist = 1;
        output_check_script<char, unsigned> ob(ASVECTOR(seq1), ASVECTOR(seq2));
        unsigned d = edit_distance(seq1, seq2, _script=ob, _substitution=true_type());
        ob.finalize(dist);
        BOOST_CHECK_MESSAGE(ob.correct, "incorrect edit script: '" << ob.ss.str() << "'  seq1='" << seq1 << "'  seq2='" << seq2 << "'");
        BOOST_CHECK_MESSAGE(d == dist, "incorrect edit distance " << d << "(expected " << dist << ")  seq1='" << seq1 << "' seq2='" << seq2 << "'  script='" << ob.ss.str() <<"'");
    }
    {
        const std::string seq1 = "abc";
        const std::string seq2 = "aBc";
        const unsigned dist = 0;
        output_check_script<char, unsigned, case_equal> ob(ASVECTOR(seq1), ASVECTOR(seq2));
        unsigned d = edit_distance(seq1, seq2, _script=ob, _substitution=true_type(), _equal=case_equal());
        ob.finalize(dist);
        BOOST_CHECK_MESSAGE(ob.correct, "incorrect edit script: '" << ob.ss.str() << "'  seq1='" << seq1 << "'  seq2='" << seq2 << "'");
        BOOST_CHECK_MESSAGE(d == dist, "incorrect edit distance " << d << "(expected " << dist << ")  seq1='" << seq1 << "' seq2='" << seq2 << "'  script='" << ob.ss.str() <<"'");
    }
}

BOOST_AUTO_TEST_CASE(myers_empty) {
    CHECK_EDIT_ALIGNMENT("", "", 0);
}

BOOST_AUTO_TEST_CASE(myers_equal) {
    CHECK_EDIT_ALIGNMENT("a", "a", 0);
    CHECK_EDIT_ALIGNMENT("aa", "aa", 0);
    CHECK_EDIT_ALIGNMENT("aaa", "aaa", 0);
}

BOOST_AUTO_TEST_CASE(myers_basis) {
    CHECK_EDIT_ALIGNMENT("", "a", 1);
    CHECK_EDIT_ALIGNMENT("", "aa", 2);
    CHECK_EDIT_ALIGNMENT("", "aaa", 3);

    CHECK_EDIT_ALIGNMENT("a", "", 1);
    CHECK_EDIT_ALIGNMENT("aa", "", 2);
    CHECK_EDIT_ALIGNMENT("aaa", "", 3);
}

BOOST_AUTO_TEST_CASE(myers_prefix_suffix) {
    CHECK_EDIT_ALIGNMENT("a", "aa", 1);
    CHECK_EDIT_ALIGNMENT("a", "aaa", 2);
    CHECK_EDIT_ALIGNMENT("a", "aaaa", 3);
    CHECK_EDIT_ALIGNMENT("a", "aaaaa", 4);

    CHECK_EDIT_ALIGNMENT("aa", "a", 1);
    CHECK_EDIT_ALIGNMENT("aaa", "a", 2);
    CHECK_EDIT_ALIGNMENT("aaaa", "a", 3);
    CHECK_EDIT_ALIGNMENT("aaaaa", "a", 4);
}

BOOST_AUTO_TEST_CASE(myers_no_equal) {
    CHECK_EDIT_ALIGNMENT("a", "x", 2);
    CHECK_EDIT_ALIGNMENT("ab", "xy", 4);
    CHECK_EDIT_ALIGNMENT("abc", "xyz", 6);

    CHECK_EDIT_ALIGNMENT("a", "xy", 3);
    CHECK_EDIT_ALIGNMENT("a", "xyz", 4);

    CHECK_EDIT_ALIGNMENT("ab", "x", 3);
    CHECK_EDIT_ALIGNMENT("abc", "x", 4);
}

BOOST_AUTO_TEST_CASE(myers_equal_runs) {
    CHECK_EDIT_ALIGNMENT("aqc", "arc", 2);
    CHECK_EDIT_ALIGNMENT("aqc", "xqz", 4);

    CHECK_EDIT_ALIGNMENT("aqqc", "arrc", 4);
    CHECK_EDIT_ALIGNMENT("aqqc", "xqqz", 4);

    CHECK_EDIT_ALIGNMENT("ax", "abx", 1);
    CHECK_EDIT_ALIGNMENT("abx", "ax", 1);

    CHECK_EDIT_ALIGNMENT("ax", "abbx", 2);
    CHECK_EDIT_ALIGNMENT("abbx", "ax", 2);
}

BOOST_AUTO_TEST_CASE(interior_0) {
    CHECK_EDIT_ALIGNMENT("axc", "mxn", 4u);
    CHECK_EDIT_ALIGNMENT("abxc", "mnxo", 6u);
    CHECK_EDIT_ALIGNMENT("axbc", "mxno", 6u);
    CHECK_EDIT_ALIGNMENT("abxc", "mxno", 6u);
    CHECK_EDIT_ALIGNMENT("axbc", "mnxo", 6u);
}

BOOST_AUTO_TEST_CASE(interior_1) {
    CHECK_EDIT_ALIGNMENT("a", "b", 2u);
    CHECK_EDIT_ALIGNMENT("ab", "cd", 4u);
    CHECK_EDIT_ALIGNMENT("abc", "def", 6u);
    CHECK_EDIT_ALIGNMENT("abcd", "efgh", 8u);
    CHECK_EDIT_ALIGNMENT("abcde", "fghij", 10u);
    CHECK_EDIT_ALIGNMENT("abcdef", "ghijkl", 12u);
}

BOOST_AUTO_TEST_CASE(interior_2) {
    CHECK_EDIT_ALIGNMENT("a", "bc", 3u);
    CHECK_EDIT_ALIGNMENT("a", "bcd", 4u);
    CHECK_EDIT_ALIGNMENT("a", "bcde", 5u);
    CHECK_EDIT_ALIGNMENT("a", "bcdef", 6u);

    CHECK_EDIT_ALIGNMENT("ab", "cde", 5u);
    CHECK_EDIT_ALIGNMENT("ab", "cdef", 6u);
    CHECK_EDIT_ALIGNMENT("ab", "cdefg", 7u);
    CHECK_EDIT_ALIGNMENT("ab", "cdefgh", 8u);

    CHECK_EDIT_ALIGNMENT("abc", "defg", 7u);
    CHECK_EDIT_ALIGNMENT("abc", "defgh", 8u);
    CHECK_EDIT_ALIGNMENT("abc", "defghi", 9u);
    CHECK_EDIT_ALIGNMENT("abc", "defghij", 10u);
}

BOOST_AUTO_TEST_CASE(interior_2_switch) {
    CHECK_EDIT_ALIGNMENT("ab", "m", 3u);
    CHECK_EDIT_ALIGNMENT("abc", "m", 4u);
    CHECK_EDIT_ALIGNMENT("abcd", "m", 5u);
    CHECK_EDIT_ALIGNMENT("abcde", "m", 6u);

    CHECK_EDIT_ALIGNMENT("ab", "mn", 4u);
    CHECK_EDIT_ALIGNMENT("abc", "mn", 5u);
    CHECK_EDIT_ALIGNMENT("abcd", "mn", 6u);
    CHECK_EDIT_ALIGNMENT("abcde", "mn", 7u);

    CHECK_EDIT_ALIGNMENT("ab", "mno", 5u);
    CHECK_EDIT_ALIGNMENT("abc", "mno", 6u);
    CHECK_EDIT_ALIGNMENT("abcd", "mno", 7u);
    CHECK_EDIT_ALIGNMENT("abcde", "mno", 8u);
}

BOOST_AUTO_TEST_CASE(interior_3) {
    CHECK_EDIT_ALIGNMENT("axbyc", "mxnyo", 6u);
    CHECK_EDIT_ALIGNMENT("axbyc", "mxnyop", 7u);
    CHECK_EDIT_ALIGNMENT("axbyc", "mxnyopq", 8u);
    CHECK_EDIT_ALIGNMENT("axbyc", "mxnyopqr", 9u);

    CHECK_EDIT_ALIGNMENT("axbyzc", "mxnyzo", 6u);
    CHECK_EDIT_ALIGNMENT("axbyzc", "mxnyzop", 7u);
    CHECK_EDIT_ALIGNMENT("axbyzc", "mxnyzopq", 8u);
    CHECK_EDIT_ALIGNMENT("axbyzc", "mxnyzopqr", 9u);

    CHECK_EDIT_ALIGNMENT("awxbyzc", "mwxnyzo", 6u);
    CHECK_EDIT_ALIGNMENT("awxbyzc", "mwxnyzop", 7u);
    CHECK_EDIT_ALIGNMENT("awxbyzc", "mwxnyzopq", 8u);
    CHECK_EDIT_ALIGNMENT("awxbyzc", "mwxnyzopqr", 9u);
}

BOOST_AUTO_TEST_CASE(interior_3_switch) {
    CHECK_EDIT_ALIGNMENT("axbyc", "mxnyo", 6u);
    CHECK_EDIT_ALIGNMENT("axbycd", "mxnyo", 7u);
    CHECK_EDIT_ALIGNMENT("axbycde", "mxnyo", 8u);
    CHECK_EDIT_ALIGNMENT("axbycdef", "mxnyo", 9u);

    CHECK_EDIT_ALIGNMENT("axbyzc", "mxnyzo", 6u);
    CHECK_EDIT_ALIGNMENT("axbyzcd", "mxnyzo", 7u);
    CHECK_EDIT_ALIGNMENT("axbyzcde", "mxnyzo", 8u);
    CHECK_EDIT_ALIGNMENT("axbyzcdef", "mxnyzo", 9u);

    CHECK_EDIT_ALIGNMENT("awxbyzc", "mwxnyzo", 6u);
    CHECK_EDIT_ALIGNMENT("awxbyzcd", "mwxnyzo", 7u);
    CHECK_EDIT_ALIGNMENT("awxbyzcde", "mwxnyzo", 8u);
    CHECK_EDIT_ALIGNMENT("awxbyzcdef", "mwxnyzo", 9u);
}


BOOST_AUTO_TEST_CASE(interior_4) {
    CHECK_EDIT_ALIGNMENT("waaa1bbbx", "yaaabbbzzz", 7u);
    CHECK_EDIT_ALIGNMENT("waaa12bbbx", "yaaabbbzzz", 8u);
    CHECK_EDIT_ALIGNMENT("waaa123bbbx", "yaaabbbzzzz", 10u);
}

BOOST_AUTO_TEST_CASE(interior_4_switch) {
    CHECK_EDIT_ALIGNMENT("waaa1bbbxxx", "yaaabbbz", 7u);
    CHECK_EDIT_ALIGNMENT("waaa12bbbxxx", "yaaabbbz", 8u);
}

BOOST_AUTO_TEST_CASE(interior_5) {
    CHECK_EDIT_ALIGNMENT("xxa", "mxx", 2u);
    CHECK_EDIT_ALIGNMENT("xxab", "mnxx", 4u);
    CHECK_EDIT_ALIGNMENT("xxabc", "mnoxx", 6u);
    CHECK_EDIT_ALIGNMENT("xxabcd", "mnopxx", 8u);
    CHECK_EDIT_ALIGNMENT("xxabc", "mnopxx", 7u);
}

BOOST_AUTO_TEST_CASE(interior_6) {
    CHECK_EDIT_ALIGNMENT("mxxxxn", "axxxbxx", 5u);
    CHECK_EDIT_ALIGNMENT("axxxbxx", "mxxxxn", 5u);
}

BOOST_AUTO_TEST_CASE(crosscheck_error_0) {
    CHECK_EDIT_ALIGNMENT("b2", "xx8a4b", 6u);
}

BOOST_AUTO_TEST_CASE(crosscheck_error_1) {
//    string seq1= "xxxxxxxxxxxxxxx8a4bd5gj9i87ddxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx4882ijxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";
//    string seq2= "xxxxxxxxxxxxxb231cjgd3hi0c5g887d9exxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";

//    string seq1= "xxxx8a4bd5gj9i87ddxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx4882ijxx";
//    string seq2= "xxb231cjgd3hi0c5g887d9exxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";

    string seq1= "xxxx8a4bd5gj9i87ddxxxxxxxx";
    string seq2= "xxb231cjgd3hi0c5g887d9exxx";

    unsigned d = edit_distance(seq1, seq2);
    std::cout << "dist= " << d << std::endl;

    CHECK_EDIT_ALIGNMENT(seq1, seq2, d);
}

BOOST_AUTO_TEST_CASE(allow_sub_1) {
    CHECK_EDIT_ALIGNMENT_ARG("abc", "xyz", _substitution=true, 3);
    CHECK_EDIT_ALIGNMENT_ARG("abc", "xyz", _substitution=false, 6);

    CHECK_EDIT_ALIGNMENT_ARG("abc", "xyz", _substitution=boost::true_type(), 3);
    CHECK_EDIT_ALIGNMENT_ARG("abc", "xyz", _substitution=boost::false_type(), 6);

    CHECK_EDIT_ALIGNMENT_ARG("aqc", "xqz", _substitution=boost::true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG("aqc", "xqz", _substitution=boost::false_type(), 4);

    CHECK_EDIT_ALIGNMENT_ARG("aqcr", "xqzr", _substitution=boost::true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG("aqcr", "xqzr", _substitution=boost::false_type(), 4);

    CHECK_EDIT_ALIGNMENT_ARG("raqc", "rxqz", _substitution=boost::true_type(), 2);
    CHECK_EDIT_ALIGNMENT_ARG("raqc", "rxqz", _substitution=boost::false_type(), 4);
}

#if 0
BOOST_AUTO_TEST_CASE(max_cost_0) {
    CHECK_EDIT_ALIGNMENT_2ARG("abc", "xyz", _substitution=true_type(), _max_cost=0, 3);
    CHECK_EDIT_ALIGNMENT_2ARG("abc", "xyz", _substitution=true_type(), _max_cost=1, 3);
    CHECK_EDIT_ALIGNMENT_2ARG("abc", "xyz", _substitution=true_type(), _max_cost=2, 3);
    CHECK_EDIT_ALIGNMENT_2ARG("abc", "xyz", _substitution=true_type(), _max_cost=3, 3);
    CHECK_EDIT_ALIGNMENT_2ARG("abc", "xyz", _substitution=true_type(), _max_cost=4, 3);

    CHECK_EDIT_ALIGNMENT_2ARG("abc", "xyz", _cost=unit_cost_test(), _max_cost=0, 6);
    CHECK_EDIT_ALIGNMENT_2ARG("abc", "xyz", _cost=unit_cost_test(), _max_cost=1, 6);
    CHECK_EDIT_ALIGNMENT_2ARG("abc", "xyz", _cost=unit_cost_test(), _max_cost=2, 6);
    CHECK_EDIT_ALIGNMENT_2ARG("abc", "xyz", _cost=unit_cost_test(), _max_cost=3, 6);
    CHECK_EDIT_ALIGNMENT_2ARG("abc", "xyz", _cost=unit_cost_test(), _max_cost=4, 6);

    CHECK_EDIT_ALIGNMENT_ARG("abc", "xyz", _max_cost=0, 6);
    CHECK_EDIT_ALIGNMENT_ARG("abc", "xyz", _max_cost=1, 6);
    CHECK_EDIT_ALIGNMENT_ARG("abc", "xyz", _max_cost=2, 6);
    CHECK_EDIT_ALIGNMENT_ARG("abc", "xyz", _max_cost=3, 6);
    CHECK_EDIT_ALIGNMENT_ARG("abc", "xyz", _max_cost=4, 6);
}

BOOST_AUTO_TEST_CASE(max_cost_1) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 100;
    random_localized_deviations(seqdata, N, 100, 5, 10);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            output_check_script_long_string out(seqdata[i], seqdata[j]);

            unsigned int d = edit_distance(seqdata[i], seqdata[j], _script=out);
            if (d < 2) continue;

            unsigned int dub = seqdata[i].size() + seqdata[j].size();
            unsigned int dt;

            out.reset();
            dt = edit_distance(seqdata[i], seqdata[j], _script=out, _max_cost=d/2);
            out.finalize(dt);
            BOOST_CHECK_GT(dt, d/2);
            BOOST_CHECK_GE(dt, d);
            BOOST_CHECK_LE(dt, dub);
            BOOST_CHECK_MESSAGE(out.correct, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j] <<"'\n\n");

            out.reset();
            dt = edit_distance(seqdata[i], seqdata[j], _script=out, _max_cost=d-1);
            out.finalize(dt);
            BOOST_CHECK_GT(dt, d-1);
            BOOST_CHECK_GE(dt, d);
            BOOST_CHECK_LE(dt, dub);
            BOOST_CHECK_MESSAGE(out.correct, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j] <<"'\n\n");

            out.reset();
            dt = edit_distance(seqdata[i], seqdata[j], _script=out, _max_cost=d);
            out.finalize(dt);
            BOOST_CHECK_EQUAL(dt, d);
            BOOST_CHECK(out.correct);

            out.reset();
            dt = edit_distance(seqdata[i], seqdata[j], _script=out, _max_cost=d+1);
            out.finalize(dt);
            BOOST_CHECK_EQUAL(dt, d);
            BOOST_CHECK(out.correct);

            out.reset();
            dt = edit_distance(seqdata[i], seqdata[j], _script=out, _max_cost=d+2);
            out.finalize(dt);
            BOOST_CHECK_EQUAL(dt, d);
            BOOST_CHECK(out.correct);

            out.reset();
            BOOST_CHECK_THROW(edit_distance(seqdata[i], seqdata[j], _script=out, _max_cost=d-1, _max_cost_exception=true), max_edit_cost_exception);

            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n");
}


BOOST_AUTO_TEST_CASE(max_cost_error_1) {
    std::string seq1 = "xxxxxjja0xxxxxxxxxxxxb544b0xxxxxxxxxxxxxxxxx2xxxxxxxxxxxxxxxxxxxxxxxicff76xxxxxxxxxxxxxxxdxxxxxxxxxx";
    std::string seq2 = "xxxxxxx16g7xxxxxxxxxxxxxxx8f5c5xxxxxxxxxxxxxah9xxxxxxxxxxxxxxbxxxxxxxxxxxxxxxxxxxxxxxxxxxj5432539cxxxxxxx";
    output_check_script_long_string out(seq1, seq2);

    int d = edit_distance(seq1, seq2, _script=out, _max_cost=20);
    out.finalize(d);
    BOOST_CHECK_MESSAGE(out.correct, "\n\nseq1= '" << seq1 << "'\nseq2= '"<< seq2 <<"'\n\n");
}


BOOST_AUTO_TEST_CASE(max_cost_2) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 100;
    random_localized_deviations(seqdata, N, 100, 5, 10);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            output_check_script_long_string out(seqdata[i], seqdata[j]);

            unsigned int d = edit_distance(seqdata[i], seqdata[j], _script=out, _substitution=true_type());
            if (d < 2) continue;

            unsigned int dub = std::max(seqdata[i].size(), seqdata[j].size());
            unsigned int dt;

            out.reset();
            dt = edit_distance(seqdata[i], seqdata[j], _script=out, _substitution=true_type(), _max_cost=d/2);
            out.finalize(dt);
            BOOST_CHECK_GT(dt, d/2);
            BOOST_CHECK_GE(dt, d);
            BOOST_CHECK_LE(dt, dub);
            BOOST_CHECK(out.correct);

            out.reset();
            dt = edit_distance(seqdata[i], seqdata[j], _script=out, _substitution=true_type(), _max_cost=d-1);
            out.finalize(dt);
            BOOST_CHECK_GT(dt, d-1);
            BOOST_CHECK_GE(dt, d);
            BOOST_CHECK_LE(dt, dub);
            BOOST_CHECK(out.correct);

            out.reset();
            dt = edit_distance(seqdata[i], seqdata[j], _script=out, _substitution=true_type(), _max_cost=d);
            out.finalize(dt);
            BOOST_CHECK_EQUAL(dt, d);
            BOOST_CHECK(out.correct);

            out.reset();
            dt = edit_distance(seqdata[i], seqdata[j], _script=out, _substitution=true_type(), _max_cost=d+1);
            out.finalize(dt);
            BOOST_CHECK_EQUAL(dt, d);
            BOOST_CHECK(out.correct);

            out.reset();
            dt = edit_distance(seqdata[i], seqdata[j], _script=out, _substitution=true_type(), _max_cost=d+2);
            out.finalize(dt);
            BOOST_CHECK_EQUAL(dt, d);
            BOOST_CHECK(out.correct);

            out.reset();
            BOOST_CHECK_THROW(edit_distance(seqdata[i], seqdata[j], _script=out, _substitution=true_type(), _max_cost=d-1, _max_cost_exception=true), max_edit_cost_exception);

            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n");
}

#endif


BOOST_AUTO_TEST_CASE(long_sequences) {
    CHECK_EDIT_ALIGNMENT_ARG("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", 
                         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", _substitution=true_type(), 
                         0);
    CHECK_EDIT_ALIGNMENT_ARG("xxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", 
                         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxx", _substitution=true_type(), 
                         2);
    CHECK_EDIT_ALIGNMENT_ARG("xxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", 
                         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*xxxxxxxxxxxxxxxxxxxxxxxxxx", _substitution=true_type(), 
                         4);
    CHECK_EDIT_ALIGNMENT_ARG("xxxxxxxxxxxxxxxxxxxxxxxx**xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx**xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", 
                         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx**xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx**xxxxxxxxxxxxxxxxxxxxxxxxx", _substitution=true_type(), 
                         8);
}


BOOST_AUTO_TEST_CASE(failure_1) {
    std::string seq1 = "xxxxxxxxx3d07a05d385h77xxxxxxxxxxxx";
    std::string seq2 = "xbd9a3d2gjf6b7a77hjcxxxxxxxxxxxxxxx";
    output_check_script_string out(seq1, seq2);
    edit_distance(seq1, seq2, _script = out, _substitution=true_type(), _cost = cost_mixed_ops());
    out.finalize();
    BOOST_CHECK_MESSAGE(out.correct, "\n\nseq1= '" << seq1 << "'\nseq2= '"<< seq2 <<"'\n\n");
}



BOOST_AUTO_TEST_CASE(timing_1) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 100;
    random_localized_deviations(seqdata, N, 100000, 5, 20, 100);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            output_check_script_long_string out(seqdata[i], seqdata[j]);
            unsigned int d = edit_distance(seqdata[i], seqdata[j], _script = out, _substitution=true_type(), _cost = cost_mixed_ops());
            out.finalize(d);
            BOOST_CHECK(out.correct);
            BOOST_CHECK(d <= 2*seqdata[i].size());
            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n" );
}

BOOST_AUTO_TEST_CASE(crosscheck_1) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 1000;
    random_localized_deviations(seqdata, N, 100, 2, 25);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            output_check_script_long_string out(seqdata[i], seqdata[j]);
            unsigned int d1 = edit_distance(seqdata[i], seqdata[j], _script = out, _substitution=true_type(), _cost = cost_mixed_ops());
            unsigned int d2 = edit_distance(seqdata[i], seqdata[j], _substitution=true_type(), _cost = cost_mixed_ops());
            out.finalize(d2);
            // verify that the edit script is a correct one: it transforms seq1 into seq2
            BOOST_CHECK_MESSAGE(out.correct, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j] <<"'\n\n");
            // cross-check 
            BOOST_CHECK_MESSAGE(d1==d2, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j] <<"'\n\n");
            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n" );
}

BOOST_AUTO_TEST_CASE(myers_sssp_crosscheck_1) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 100;
    random_localized_deviations(seqdata, N, 100, 2, 25);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            output_check_script_long_string out1(seqdata[i], seqdata[j]);
            output_check_script_long_string out2(seqdata[i], seqdata[j]);
            // Myers algorithm
            unsigned int d1 = edit_distance(seqdata[i], seqdata[j], _script=out1, _substitution=boost::false_type());
            out1.finalize(d1);
            // SSSP algorithm
            unsigned int d2 = edit_distance(seqdata[i], seqdata[j], _script=out2, _substitution=boost::false_type(), _cost = unit_cost_test());
            out2.finalize(d2);
            // verify that the edit script is a correct one: it transforms seq1 into seq2
            BOOST_CHECK_MESSAGE(out1.correct, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j] <<"'\n\n");
            BOOST_CHECK_MESSAGE(out2.correct, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j] <<"'\n\n");
            // cross-check 
            BOOST_CHECK_MESSAGE(d1==d2, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j] <<"'\n\n");
            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n" );
}

BOOST_AUTO_TEST_CASE(myers_dist_path_crosscheck_1) {
    srand(time(0));
    vector<std::string> seqdata;
    const int N = 100000;
    random_localized_deviations(seqdata, N, 100, 5, 10);
    int n = 0;
    double t0 = time(0);
    for (unsigned int i = 0;  i < seqdata.size();  ++i) {
        if (n >= N) break;
        for (unsigned int j = 0;  j < i;  ++j) {
            BOOST_TEST_CHECKPOINT("n= " << n << "   i= " << i << "   j= " << j);
            output_check_script_long_string out(seqdata[i], seqdata[j]);
            unsigned int d2 = edit_distance(seqdata[i], seqdata[j], _substitution=boost::false_type());
                //BOOST_TEST_CHECKPOINT("aaa");
            unsigned int d1 = edit_distance(seqdata[i], seqdata[j], _script = out, _substitution=boost::false_type());
                //BOOST_TEST_CHECKPOINT("bbb");
            out.finalize(d2);
            // verify that the edit script is a correct one: it transforms seq1 into seq2
            BOOST_CHECK_MESSAGE(out.correct, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j] <<"'\n\n");
            // cross-check 
            BOOST_CHECK_MESSAGE(d1==d2, "\n\nseq1= '" << seqdata[i] << "'\nseq2= '"<< seqdata[j] <<"'\n\n");
            if (++n >= N) break;
        }
    }
    double tt = time(0) - t0;
    BOOST_TEST_MESSAGE("time= " << tt << " sec   n= " << n << "   mean-time= " << tt/double(n) << "\n");
}

BOOST_AUTO_TEST_SUITE_END()
