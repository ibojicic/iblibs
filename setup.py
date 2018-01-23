from setuptools import setup, find_packages

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='iblibs',
    version='1.1.0',
    description='Function Libraries',
    long_description=readme,
    author='Ivan Bojicic',
    author_email='ibojicic@gmail.com',
    license=license,
    packages=find_packages(),
    install_requires=[
        # 'Click',
    ],
    entry_points={
        'console_scripts': [
            'python_driver=qb_drivers.python_driver:cli',
        ]
    }

)
