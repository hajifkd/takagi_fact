import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="takagi_fact",
    version="0.2.1",
    author="Hajime Fukuda",
    author_email="hajime.fukuda@me.com",
    description="A libary for the symmetric SVD and the Takagi factorization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hajifkd/takagi_fact/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    install_requires=[
        'mpmath',
    ]
)
