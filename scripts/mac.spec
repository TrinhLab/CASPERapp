# -*- mode: python ; coding: utf-8 -*-
import sys
sys.setrecursionlimit(5000)
block_cipher = None


a = Analysis(['main.py'],
             pathex=['/Users/ddooley/bioinformatics_packages/individual_packages/crispr_tools/CASPERapp'],
             datas=[
				('OffTargetFolder', 'OffTargetFolder'),
				('SeqFinderFolder', 'SeqFinderFolder'),
				('*.ui','.'),
				('CASPERinfo', '.'),
				('CASPER-logo.jpg', '.'),
				('CASPER_icon.icns', '.'),
				('genomeBrowserTemplate.html', '.'),
                ('/Users/ddooley/venvs/casper_env/lib/python3.8/site-packages/azimuth','azimuth')
			 ],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[
				 'cryptography',
				 'Crypto',
                 'tkinter'
			 ],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='CASPERapp',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False)
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='CASPERapp')

app = BUNDLE(coll, icon='CASPER_icon.icns',
             name='CASPERapp.app',
             version='2.0.1',
             bundle_identifier=None)
