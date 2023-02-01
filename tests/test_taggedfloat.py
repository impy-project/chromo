from chromo.util import TaggedFloat


def test_taggedfloat():
    tf = TaggedFloat(1)
    assert isinstance(tf, TaggedFloat)
    assert not isinstance(tf, float)
    assert repr(tf) == "TaggedFloat(1.0)"
    assert tf == TaggedFloat(1)
    assert tf * 2 == TaggedFloat(2)
    assert 2 * tf == TaggedFloat(2)
    assert tf + 1 == TaggedFloat(2)
    assert 1 + tf == TaggedFloat(2)
    assert tf - 1 == TaggedFloat(0)
    assert 1 - tf == TaggedFloat(0)
    assert tf / 2 == TaggedFloat(0.5)
    assert 2 / tf == TaggedFloat(2)
    assert float(tf) == 1


class Foo(TaggedFloat):
    pass


def test_derived():
    d = Foo(1)
    assert isinstance(d, Foo)
    assert not isinstance(d, float)
    assert repr(d) == "Foo(1.0)"
    assert d == Foo(1)
    assert d * 2 == Foo(2)
    assert 2 * d == Foo(2)
    assert d + 1 == Foo(2)
    assert 1 + d == Foo(2)
    assert d - 1 == Foo(0)
    assert 1 - d == Foo(0)
    assert d / 2 == Foo(0.5)
    assert 2 / d == Foo(2)
    assert float(d) == 1
