from skbuild import setup

setup(
    packages=['overlap'],
    package_dir={'': 'python'},
    cmake_install_dir='python/overlap',
)
