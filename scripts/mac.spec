block_cipher = None

a = Analysis(['src/main.py'],
             pathex=['src'],
             datas=[
                ('assets', 'assets'),
                ('config', 'config'),
                ('logs', 'logs'),
                ('src', 'src'),
                ('genomeBrowserTemplate.html', '.'),
             ],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
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
          console=False,
          disable_windowed_traceback=False,
          target_arch=None,
          codesign_identity=None,
          entitlements_file=None,
          icon='assets/CASPER_icon.icns')
          
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='CASPERapp')

app = BUNDLE(coll, 
             name='CASPERapp.app',
             icon='assets/CASPER_icon.icns',
             version='2.0.1',
             bundle_identifier=None)

# 1. Have the mac.spec in the app directory
# 2. pyinstaller mac.spec
# 3. mkdir -p dist/dmg
# 4. rm -r dist/dmg/*
# 5. Manual copy of the app into dist/dmg
# 6. create-dmg \
#   --volname "CASPERapp" \
#   --window-pos 200 120 \
#   --window-size 600 300 \
#   --icon-size 100 \
#   --icon "CASPERapp.app" 175 120 \
#   --hide-extension "CASPERapp.app" \
#   --app-drop-link 425 120 \
#   "dist/CASPERapp.dmg" \
#   "dist/dmg/"