# -*- mode: python ; coding: utf-8 -*-
import sys
sys.setrecursionlimit(5000)
block_cipher = None


a = Analysis(['main.py'],
             pathex=['C:/Users/Trinh^ Lab/Desktop/CASPERapp'],
             datas=[
				('OffTargetFolder', 'OffTargetFolder'),
				('SeqFinderFolder', 'SeqFinderFolder'),
				('*.ui','.'),
				('CASPERinfo', '.'),
				('CASPER-logo.jpg', '.'),
				('CASPER_icon.icns', '.'),
				('CASPER_icon.ico', '.'),
				('genomeBrowserTemplate.html', '.'),
				('C:/Users/Trinh Lab/Desktop/casper_env/Lib/site-packages/azimuth', 'azimuth')
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
