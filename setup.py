from setuptools import setup, find_packages

setup(
    name='arca',
    version='0.1',
    author='Denis Fitz',
    description='A pattern language for music',
    packages=find_packages(),
    install_requires=[
        'random',
        'copy',
        'DTXML',
    ],
)
