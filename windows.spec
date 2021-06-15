# -*- mode: python ; coding: utf-8 -*-
import sys
sys.setrecursionlimit(5000)
block_cipher = None


a = Analysis(['main.py'],
             pathex=['C:/Users/Tfry/Desktop/CASPERapp'],
             datas=[
				('OffTargetFolder', 'OffTargetFolder'),
				('SeqFinderFolder', 'SeqFinderFolder'),
				('CASPERinfo', '.'),
            ('settings_icon.png', '.'),
				('Casperimg.icns', '.'),
				('Codonimg.icns', '.'),
				('myicon.ico', '.'),
				('mainart.jpg', '.'),
				('cas9image.ico', '.'),
				('*.ui', '.'),
				('CASPER-logo.jpg', '.'),
				('genomeBrowserTemplate.html', '.')
			 ],
             hiddenimports=[
				'repeats_vs_seeds_line', 
				'seeds_vs_repeats_bar', 
				'repeats_vs_chromo',
				'pop_analysis_3dgraph',
				'pop_analysis_repeats_graph',
				'pop_analysis_venn_diagram',
				'pkg_resources.py2_warn'
			],
             hookspath=[],
             runtime_hooks=[],
             excludes=[
				 'cryptography',
				 'Crypto'
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
