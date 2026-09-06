"""Unit tests for the queries.txt parser used by scripts/verify_counts.py.

These are pure-text tests: no Hasura, no network, no database.
"""
import pytest

from scripts.verify_counts import Query, extract_counts, parse_queries


class TestParseQueries:
    def test_single_query(self):
        text = """
query Simple {
  Simple: variants_aggregate { aggregate { count } }
}
"""
        (q,) = parse_queries(text)
        assert q.name == "Simple"
        assert q.text.startswith("query Simple {")
        assert q.text.rstrip().endswith("}")
        assert q.has_variables is False
        assert q.expect == {}

    def test_multiple_queries_are_split(self):
        text = """
query A {
  A: variants_aggregate { aggregate { count } }
}

query B {
  B: genes_aggregate { aggregate { count } }
}
"""
        names = [q.name for q in parse_queries(text)]
        assert names == ["A", "B"]

    def test_nested_braces_do_not_end_the_block_early(self):
        text = """
query Deep {
  Deep: variants_aggregate(where: {a: {b: {c: {_eq: "x"}}}}) { aggregate { count } }
}

query After {
  After: genes_aggregate { aggregate { count } }
}
"""
        queries = parse_queries(text)
        assert [q.name for q in queries] == ["Deep", "After"]
        assert "After" not in queries[0].text

    def test_braces_inside_comments_are_ignored(self):
        """A commented-out filter line contains an unbalanced brace count; the block
        must still terminate at the right place. This is exactly the shape of the
        previously commented-out off-target line in queries.txt."""
        text = """
query Commented {
  # variants_guides: {edit_type: {_eq: "Rescue"}, guide: {hits_85: {_eq: 0}}},
  Commented: variants_aggregate { aggregate { count } }
}

query After {
  After: genes_aggregate { aggregate { count } }
}
"""
        queries = parse_queries(text)
        assert [q.name for q in queries] == ["Commented", "After"]

    def test_braces_inside_strings_are_ignored(self):
        text = """
query Stringy {
  Stringy: variants_aggregate(where: {name: {_eq: "a{b"}}) { aggregate { count } }
}

query After {
  After: genes_aggregate { aggregate { count } }
}
"""
        assert [q.name for q in parse_queries(text)] == ["Stringy", "After"]

    def test_query_with_variables_is_flagged(self):
        text = """
query Paged($limit: Int!, $offset: Int!) {
  Paged: variants(limit: $limit, offset: $offset) { id }
}
"""
        (q,) = parse_queries(text)
        assert q.name == "Paged"
        assert q.has_variables is True

    def test_no_space_before_brace(self):
        """queries.txt contains `query NonG2A_StopGained_Pass{` with no space."""
        text = """
query Tight{
  Tight: variants_aggregate { aggregate { count } }
}
"""
        (q,) = parse_queries(text)
        assert q.name == "Tight"

    def test_line_number_is_recorded(self):
        text = "\n\nquery A {\n  A: x_aggregate { aggregate { count } }\n}\n"
        (q,) = parse_queries(text)
        assert q.line == 3


class TestExpectAnnotations:
    def test_bare_expect_binds_to_the_single_alias(self):
        text = """
query Solo {
  # expect: 2421
  Solo: variants_aggregate { aggregate { count } }
}
"""
        (q,) = parse_queries(text)
        assert q.expect == {None: 2421}

    def test_aliased_expect(self):
        text = """
query Funnel {
  # expect TotalVariants: 329279
  # expect SNVsVariants: 281812
  TotalVariants: variants_aggregate { aggregate { count } }
  SNVsVariants: variants_aggregate { aggregate { count } }
}
"""
        (q,) = parse_queries(text)
        assert q.expect == {"TotalVariants": 329279, "SNVsVariants": 281812}

    def test_expect_tolerates_thousands_separators_and_spacing(self):
        text = """
query Solo {
  #expect:  2,421
  Solo: variants_aggregate { aggregate { count } }
}
"""
        (q,) = parse_queries(text)
        assert q.expect == {None: 2421}

    def test_no_annotation_yields_empty_expect(self):
        text = "query Solo {\n  Solo: variants_aggregate { aggregate { count } }\n}\n"
        (q,) = parse_queries(text)
        assert q.expect == {}


class TestExtractCounts:
    def test_flat_aggregate(self):
        data = {"Solo": {"aggregate": {"count": 2421}}}
        assert extract_counts(data) == {"Solo": 2421}

    def test_multiple_aliases(self):
        data = {
            "TotalVariants": {"aggregate": {"count": 329279}},
            "SNVsVariants": {"aggregate": {"count": 281812}},
        }
        assert extract_counts(data) == {
            "TotalVariants": 329279,
            "SNVsVariants": 281812,
        }

    def test_nested_aggregate_under_a_relationship(self):
        data = {"subjects": {"variants": {"variants_aggregate": {"aggregate": {"count": 7}}}}}
        assert extract_counts(data) == {"subjects.variants.variants_aggregate": 7}

    def test_list_results_are_counted_by_length(self):
        data = {"rows": [{"id": 1}, {"id": 2}, {"id": 3}]}
        assert extract_counts(data) == {"rows": 3}

    def test_empty_data(self):
        assert extract_counts({}) == {}


class TestQueryVerification:
    def test_mismatch_is_reported(self):
        q = Query(name="Solo", text="", line=1, expect={None: 315}, has_variables=False)
        problems = q.check({"Solo": 404})
        assert len(problems) == 1
        assert "315" in problems[0] and "404" in problems[0]

    def test_match_is_silent(self):
        q = Query(name="Solo", text="", line=1, expect={None: 315}, has_variables=False)
        assert q.check({"Solo": 315}) == []

    def test_bare_expect_against_multiple_counts_is_ambiguous(self):
        q = Query(name="Multi", text="", line=1, expect={None: 5}, has_variables=False)
        problems = q.check({"a": 5, "b": 6})
        assert len(problems) == 1
        assert "ambiguous" in problems[0].lower()

    def test_aliased_expect_for_a_missing_alias(self):
        q = Query(name="Q", text="", line=1, expect={"nope": 1}, has_variables=False)
        problems = q.check({"other": 1})
        assert len(problems) == 1
        assert "nope" in problems[0]

    def test_no_expectation_never_fails(self):
        q = Query(name="Q", text="", line=1, expect={}, has_variables=False)
        assert q.check({"anything": 99}) == []
