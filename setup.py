from setuptools import setup, find_packages

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

    python_requires='>=3.6',
)
