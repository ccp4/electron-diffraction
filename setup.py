import setuptools

#with open("latex/README.md", "r") as fh:
#    long_description = fh.read()

setuptools.setup(
    name="ccp4ED",
    version="1.0.8",
    author="Tarik Ronan Drevon",
    author_email="ronandrevon@gmail.com",
    description="Electron diffraction utilities",
    long_description='', #long_description,
    long_description_content_type="",
    url="https://pypi.org/project/pyScatSpheres",
    project_urls={
        'Documentation': 'https://debloch.readthedocs.io/en/latest/',
        'Source':'https://github.com/ccp4/electron-diffraction',
    },

    # packages=['multislice','scattering','wallpp'],
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={'multislice/data':['splines.npy']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License ",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['numpy','scipy','matplotlib','colorama','pandas',
    'cbf','crystals','TarikDrevonUtils','easygui','tifffile','pickle5','bs4'],
)
