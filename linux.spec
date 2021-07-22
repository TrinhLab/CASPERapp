# -*- mode: python ; coding: utf-8 -*-
import sys
sys.setrecursionlimit(5000)
block_cipher = None


a = Analysis(['main.py'],
			pathex=['/home/ddooley/CASPERapp'],
			datas=[
			('OffTargetFolder', 'OffTargetFolder'),
			('SeqFinderFolder', 'SeqFinderFolder'),
			('*.ui','.'),
			('CASPERinfo', '.'),
			('Casperimg.icns', '.'),
			('Codonimg.icns', '.'),
			('myicon.ico', '.'),
			('mainart.jpg', '.'),
			('cas9image.ico', '.'),
			('settings_icon.png', '.'),
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
			a.binaries,
			a.zipfiles,
			a.datas,
			name='CASPERapp',
			debug=False,
			strip=False,
			bootloader_ignore_signals=False,
			upx=True,
			console=True)
app = BUNDLE(icon='test.icns',
			name='CASPERapp.app',
			bundle_identifier=None)
