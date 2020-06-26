import setuptools

#with open("latex/README.md", "r") as fh:
#    long_description = fh.read()

setuptools.setup(
    name="scattering",
    version="0.0.1",
    author="Tarik Ronan Drevon",
    author_email="ronandrevon@gmail.com",
    description="scattering utilities",
    long_description='', #long_description,
    long_description_content_type="",
    url="",
    packages=['scattering'],#setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License ",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    #install_requires=['matplotlib','numpy','scipy','colorama','pandas'],
)
