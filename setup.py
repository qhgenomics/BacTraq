from setuptools import setup, find_packages

setup (
    name='BacTraq',
    version='1.0',
    author="Thuy Linh Nguyen",
    author_email="ThuyLinh.Nguyen@health.qld.gov.au",
    description="A package for SNP distance clustering and consistent naming with previous cluster.",
    url= "https://github.com/qhgenomics/BacTraq",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            "bactraq=BacTraq.bactraq:main",
            "bactraq-history=BacTraq.history.makeHistory:main",
        ]
    },

    install_requires=[
        'setuptools',
        'pandas',
        'numpy',
        'scikit-learn',
        'treelib',
        'pyarrow'
    ],

    python_requires=">=3.10"
)
