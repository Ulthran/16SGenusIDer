import setuptools

setuptools.setup(
    name="GenusFinder",
    version="0.0.0",
    description="Given a bacterial 16S gene, infer the genus by placing it on a tree of similar sequences",
    author="Charlie Bushman",
    author_email="ctbushman@gmail.com",
    url="https://github.com/PennChopMicrobiomeProgram",
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    pythonpath=[".", "src"],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "idgenus=GenusFinder.command:main",
            "prepdb=GenusFinder.prepare_strain_data:main",
            "traingenus=GenusFinder.train_command:main",
        ],
    },
    install_requires=[
        "ete3",
        "eutils",
        "numpy",
        "scikit-learn",
        "tqdm",
    ],
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="GPLv2+",
)
