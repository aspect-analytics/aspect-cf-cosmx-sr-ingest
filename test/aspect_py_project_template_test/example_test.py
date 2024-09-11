import aspect_py_project_template
from aspect_py_project_template.example import get_massnumber


def test_version():
    assert aspect_py_project_template.__version__


def test_get_massnumber():
    assert get_massnumber("H2O") == 18
