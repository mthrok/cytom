from setuptools import setup

setup(
    name='cytom',
    version='0.1.0',
    packages=[
        'cytom',
    ],
    test_suite='tests',
    install_requires=[
        'GoreUtilities',
        'FlowCytometryTools',
    ],
)
