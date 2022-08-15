from setuptools import setup

with open("melt/version.py", "r") as f:
    exec(f.read())

setup(name='melt',
      version=__version__,
      description='Energy balance calculation for rocky planet interiors under the influence of tides',
      url='http://github.com/cpiaulet/melt',
      author='Caroline Piaulet',
      author_email='caroline.piaulet@umontreal.ca',
      license='GNU GPL v3.0',
      packages=['melt'],
      install_requires=['numpy', 'astropy'],
      zip_safe=False)