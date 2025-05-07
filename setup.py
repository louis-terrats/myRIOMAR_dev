from setuptools import setup, find_packages

setup(
    name='myRIOMAR',
    version='0.1.0',
    description='Satellite plume detection and analysis',
    author='Louis Terrats',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'yaml',
        # list all your deps here
    ],
    python_requires='>=3.11',
)
