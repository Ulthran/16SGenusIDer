import setuptools

setuptools.setup(
    name='16SGenusIDer',
    version="0.0.0",
    description='Given a bacterial 16S gene, infer the genus by placing it on a tree of similar sequences',
    author='Charlie Bushman',
    author_email='ctbushman@gmail.com',
    url='https://github.com/PennChopMicrobiomeProgram',
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    entry_points={
        'console_scripts': [
            'idgenus=16SGenusIDer.command:main',
        ],
    },
    install_requires=[
        'eutils',
        'vsearch',
        'dendropy',
    ],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    license='GPLv2+',
)