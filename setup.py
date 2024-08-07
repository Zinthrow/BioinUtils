import os
import codecs
import re
import sys
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

PY3 = sys.version_info[0] == 3
PY34_OR_LESS = PY3 and sys.version_info < (3, 7)


if PY34_OR_LESS:
    raise Exception("Python version less 3.7 is not supported")


def read(*parts):
    return codecs.open(os.path.join(here, *parts), 'r').read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


def _pip_requirement(req):
    if req.startswith('-r '):
        _, path = req.split()
        return reqs(*path.split('/'))
    return [req]


def _strip_comments(line):
    return line.split('#', 1)[0].strip()


def _reqs(*f):
    return [
        _pip_requirement(r) for r in (
            _strip_comments(line) for line in open('requirements.txt').readlines())
        if r]


def reqs(*f):
    """Parse requirement file.
    Example:
        reqs('requirements.txt')          # requirements.txt
    """
    return [req for subreq in _reqs(*f) for req in subreq]


def install_requires():
    """Get list of requirements required for installation."""
    return reqs('requirements.txt')


setup(
    name='bioinutils',
    version=find_version("src/bioinutils/", "__init__.py"),
    description='python package for command line blast',
    url='https://github.com/Zinthrow/BioinUtils',
    author='Alexander Larsen',
    author_email='alarsen525@gmail.com',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=install_requires(),
    python_requires=">=3.7"
)
