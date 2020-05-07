from setuptools import setup

with open('README.md') as fp:
    LONG_DESCRIPTION = fp.read()

setup(
      name='pynetstim',
      maintainer = 'Ehsan Tadyaon',
      maintainer_email = 'sunny.tadayon@gmail.com',
      version='0.2-alpha',
      description='brain network stimulation',
      long_description = LONG_DESCRIPTION,
      author='Ehsan Tadayon',
      author_email='sunny.tadayon@gmail.com',
      license='MIT',
      zip_safe=False,
      packages=['pynetstim'],
      url = 'https://github.com/EhsanTadayon/pynetstim',
      download_url='https://github.com/EhsanTadayon/pynetstim/archive/v_0.1-alpha.tar.gz',
      keywords = ['brain network','stimulation','imaging','fMRI','diffusion','TMS'],
      install_requires=[
          'numpy',
          'nipype',
          'mayavi',
          'pysurfer',
          'pandas',
          'matplotlib',
          'scipy',
          'nibabel',
          'mne',
          'nilearn',
      ],
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Developers',
          'Topic :: Scientific/Engineering',
          'Operating System :: MacOS',
          'Operating System :: Unix',
          'Topic :: Software Development :: Build Tools',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7'
          
      ],
      )


