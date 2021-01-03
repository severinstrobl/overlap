from skbuild import setup, exceptions, cmaker
from packaging.version import LegacyVersion

# Add CMake as a build requirement if no or only an outdated version is
# available.
setup_requires = []
try:
    if LegacyVersion(cmaker.get_cmake_version()) < LegacyVersion('3.10'):
        setup_requires.append('cmake')
except exceptions.SKBuildError:
    setup_requires.append('cmake')

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='overlap',
    version='0.0.5',
    author='Severin Strobl',
    author_email='git@severin-strobl.de',
    description='Exact calculation of the overlap volume and area of spheres '
                'and mesh elements',
    url='https://github.com/severinstrobl/overlap',
    license='GPLv3',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['overlap'],
    package_dir={'overlap': 'python'},
    setup_requires=setup_requires,
    install_requires=[
        'numpy'
    ],
    cmake_install_dir='python',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Topic :: Scientific/Engineering'
    ]
)
