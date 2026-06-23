; Inno Setup script for Oligolia on Windows
; Download Inno Setup from https://jrsoftware.org/isinfo.php
; Run: iscc build\inno_setup.iss

#define MyAppName "Oligolia"
#define MyAppVersion "0.3.0"
#define MyAppPublisher "Oligolia Project"
#define MyAppURL "https://github.com/moonsoup/oligolia"
#define MyAppExeName "Oligolia.exe"
#define MyAppDescription "Gene Editing and Viewing Platform"

[Setup]
AppId={{B3F7A2E1-9C4D-4E8F-A1B2-3C5D6E7F8A9B}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
AllowNoIcons=no
SourceDir=..\
OutputDir=dist
OutputBaseFilename=Oligolia-{#MyAppVersion}-Setup
SetupIconFile=assets\icon.ico
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
PrivilegesRequired=admin
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64

; Welcome / License pages
LicenseFile=LICENSE
InfoBeforeFile=build\installer_notes.txt

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked
Name: "startmenuicon"; Description: "Create Start Menu shortcut"; GroupDescription: "{cm:AdditionalIcons}"; Flags: checkedonce

[Files]
; Main application (PyInstaller output)
Source: "dist\Oligolia\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
; Start Menu
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; IconFilename: "{app}\{#MyAppExeName}"; Comment: "{#MyAppDescription}"
Name: "{group}\Uninstall {#MyAppName}"; Filename: "{uninstallexe}"
; Desktop (if user chose it)
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; IconFilename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(MyAppName, '&', '&&')}}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
Type: filesandordirs; Name: "{localappdata}\Oligolia"

[Registry]
; File associations
Root: HKCR; Subkey: ".fasta"; ValueType: string; ValueName: ""; ValueData: "OligoliaFASTA"; Flags: uninsdeletevalue
Root: HKCR; Subkey: "OligoliaFASTA"; ValueType: string; ValueName: ""; ValueData: "FASTA Sequence File"; Flags: uninsdeletekey
Root: HKCR; Subkey: "OligoliaFASTA\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\{#MyAppExeName},0"
Root: HKCR; Subkey: "OligoliaFASTA\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\{#MyAppExeName}"" ""%1"""

Root: HKCR; Subkey: ".vcf"; ValueType: string; ValueName: ""; ValueData: "OligoliaVCF"; Flags: uninsdeletevalue
Root: HKCR; Subkey: "OligoliaVCF"; ValueType: string; ValueName: ""; ValueData: "VCF Variant File"; Flags: uninsdeletekey
Root: HKCR; Subkey: "OligoliaVCF\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\{#MyAppExeName},0"
Root: HKCR; Subkey: "OligoliaVCF\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\{#MyAppExeName}"" ""%1"""
