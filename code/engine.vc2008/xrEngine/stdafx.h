#pragma once
// Modifer: sv3nk
#ifndef INGAME_EDITOR
#	define	INGAME_EDITOR
#endif // #ifndef INGAME_EDITOR

#ifndef _WIN32_WINNT
#	define _WIN32_WINNT 0x0600
#endif

#include "../xrCore/xrCore.h"
#include "../xrCore/xrAPI.h"

#ifdef _DEBUG
#	define D3D_DEBUG_INFO
#endif
#pragma warning(push)
#pragma warning(disable:4995)
#include <d3d9.h>
#pragma warning(pop)

// you must define ENGINE_BUILD then building the engine itself
// and not define it if you are about to build DLL
#ifndef NO_ENGINE_API
	#ifdef	ENGINE_BUILD
		#define DLL_API			__declspec(dllimport)
		#define ENGINE_API		__declspec(dllexport)
	#else
		#undef	DLL_API
		#define DLL_API			__declspec(dllexport)
		#define ENGINE_API		__declspec(dllimport)
	#endif
#else
	#define ENGINE_API
	#define DLL_API
#endif // NO_ENGINE_API

#define ECORE_API

// Our headers
#include "engine.h"
#include "defines.h" //TODO: replace paths defs to consts structs!
#include "consts.h"

#ifndef NO_XRLOG
#include "../xrcore/log.h"
#endif

#include "device.h"
#include "../xrcore/fs.h"
#include "../xrcdb/xrXRC.h"
#include "../xrSound/sound.h"

extern ENGINE_API CInifile *pGameIni;

#pragma comment( lib, "xrCore.lib"	)
#pragma comment( lib, "xrCDB.lib"	)
#pragma comment( lib, "xrSound.lib"	)
#pragma comment( lib, "winmm.lib"	)
#ifndef ENGINE_BUILD
#pragma comment( lib, "xrScripts.lib")
#pragma comment( lib, "d3d9.lib"	 )
#pragma comment( lib, "dinput8.lib"	 )
#pragma comment( lib, "dxguid.lib"	 )
#endif
#pragma comment( lib, "sound_static.lib")
 
#define LUABIND_DONT_COPY_STRINGS

#define READ_IF_EXISTS(ltx,method,section,name,default_value)\
	(((ltx)->line_exist(section, name)) ? ((ltx)->method(section, name)) : (default_value))
