from distutils.core import setup
import re


def get_version():
    with open('__init__.py', 'r') as f:
        content = f.read()
    p = re.compile(r'^__version__ = [\'"]([^\'\"]*)[\'"]', re.M)
    return p.search(content).group(1)


setup(name='atools',
      author='Andrew Tarzia',
      author_email='andrew.tarzia@gmail.com',
      url='https://www.github.com/andrewtarzia/atools',
      version=get_version())
