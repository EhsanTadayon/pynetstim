from setuptools import setup

setup(
      name='pynetstim',
      version='0.1',
      description='brain network stimulation',
      author='Ehsan Tadayon',
      author_email='sunny.tadayon@gmail.com',
      license='MIT',
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
          'Intended Audience :: Developers',
          'Topic :: Software Development :: Build Tools',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7'
          
      ]
      zip_safe=False)


