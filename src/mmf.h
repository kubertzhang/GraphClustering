#ifndef _MMF_H
#define _MMF_H

#include <Windows.h>
#include <WinBase.h>

#define BAD_POS 0xFFFFFFFF      // returned by SetFilePointer and GetFileSize
#define SUCCESS 0

typedef DWORD mmf_share_mode;
typedef DWORD mmf_access_mode;
typedef DWORD mmf_flags;

class MMF
{
private:
	HANDLE mmHandle;                            // 物理文件
	HANDLE mmfm;                                // 文件映射 
public:
	void mmfCreate(const char * _sharedname, const char * _filename, const DWORD _mmfsizehigh, const DWORD _mmfsizelow, float * & _mmfm_base_address);
	void mmfClose(float * & _mmfm_base_address);
};

#endif

