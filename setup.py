from setuptools import setup, find_packages, Extension
from pybind11.setup_helpers import Pybind11Extension, build_ext


# pybind11 stuff suggested by chatgpt
ext_modules = [
    Pybind11Extension(
        "fibonacci",  # Module name
        ["fibonacci.cpp"],  # Source file
    ),
]


setup(
    name='oddesa',
    version='0.0.1',
    description='Ordinary and Delay Differential Equations Solver and Analizer',
    author='Danila Bain',
    author_email='danila.bain@yandex.com',
    url='https://github.com/Kir-Dykov/oddesa',
    packages=find_packages(),
    install_requires=["regex"],                # Dependencies
    classifiers=[                       # Metadata
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    
    # pybind11 stuff suggested by chatgpt
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,

    python_requires='>=3.6',
)
