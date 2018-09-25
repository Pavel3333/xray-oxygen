#pragma once

struct GAME_DIRECTORIES {
	const char G_DATA      [12] = "$game_data$";                 //game paths
	const char G_LEVELS    [14] = "$game_levels$";
	const char G_TEXTURES  [16] = "$game_textures$";
	const char G_SHADERS   [16] = "$game_shaders$";
	const char G_SOUNDS    [14] = "$game_sounds$";
	const char G_ANIMS     [13] = "$game_anims$";
	const char G_CONFIG    [14] = "$game_config$";
	const char G_SAVES     [13] = "$game_saves$";


	const char S_ROOT      [14] = "$server_root$";               //server paths
	const char S_DATA_ROOT [19] = "$server_data_root$";


	const char IMPORT      [9]  = "$import$";                     //some other paths
	const char SOUNDS      [9]  = "$sounds$";
	const char TEXTURES    [11] = "$textures$";
	const char MAPS        [9]  = "$maps$";

	const char OBJ         [10] = "$objects$";
	const char DETAIL_OBJ  [17] = "$detail_objects$";

	const char SMOTION     [10] = "$smotion$";
	const char OMOTION     [10] = "$omotion$";
	const char OMOTIONS    [11] = "$omotions$";

	const char LEVEL       [8]  = "$level$";
	const char LEVELCACHE  [14] = "$level_cache$";


	const char LOCAL_ROOT  [13] = "$local_root$";                //roots
	const char FS_ROOT     [10] = "$fs_root$";
	const char APP_ROOT    [11] = "$app_root$";


	const char TARGET      [16] = "$target_folder$";             //some debug (???) paths
	const char BUILDCPY    [13] = "$build_copy$";


	const char TEMP        [7]  = "$temp$";                       //system paths
	const char TEMP_       [17] = "$!#%TEMP%#!$.$$$";
	const char APPDATA     [16] = "$app_data_root$";

} GameDirectories;

struct TEXTURE_FORMATS {
	const char DDS  [5] = ".dds";
	const char TGA  [5] = ".tga";
	const char THM  [5] = ".thm";
	const char BMP  [5] = ".bmp";
	const char OGM  [5] = ".ogm";
} TextureFormats;