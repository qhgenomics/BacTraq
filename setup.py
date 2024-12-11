from setuptools import setup

setup(

    scripts=['bin/bactrac'],

    entry_points={
    'console_scripts': [
        'make-history = history:makeHistory',
    ],
},
)
