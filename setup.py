from skbuild import setup, exceptions, cmaker
from packaging.version import LegacyVersion

# add CMake as a build requirement if no suitable version is available
setup_requires = []
try:
    if LegacyVersion(cmaker.get_cmake_version()) < LegacyVersion('3.12'):
        setup_requires.append('cmake')
except exceptions.SKBuildError:
    setup_requires.append('cmake')

setup(
    packages=['overlap'],
    package_dir={'overlap': 'python'},
    setup_requires=setup_requires,
    cmake_install_dir='python',
)
