#include <iostream>
#include "mmf.h"

using namespace std;


void MMF::mmfCreate(const char * _sharedname, const char * _filename, const DWORD _mmfsizehigh, const DWORD _mmfsizelow, float * & _mmfm_base_address)
{
	
	//存取模式
	mmf_access_mode access_mode = (GENERIC_READ | GENERIC_WRITE);

	//共享模式
	mmf_share_mode share_mode = 0;             // 不共享
	//mmf_share_mode share_mode = FILE_SHARE_READ | FILE_SHARE_WRITE;

	//文件属性
	mmf_flags flags = FILE_FLAG_SEQUENTIAL_SCAN;      // 针对连续访问对文件缓冲进行优化
	//mmf_flags flags = FILE_FLAG_RANDOM_ACCESS;      // 针对随机访问对文件缓冲进行优化

	DWORD error_code;

	//创建文件,返回文件句柄(HANDLE)
	mmHandle = CreateFile(_filename,		       // 要打开的文件的名字
							     access_mode,	   // GENERIC_READ 表示允许对设备进行读访问; GENERIC_WRITE 表示允许对设备进行写访问(可组合使用); 如果为零，表示只允许获取与一个设备有关的信息
		                         share_mode,	   // 0 表示不共享;FILE_SHARE_READ 和 / 或 FILE_SHARE_WRITE 表示允许对文件进行共享访问
		                         NULL,			   // SECURITY_ATTRIBUTES ，指向一个 SECURITY_ATTRIBUTES 结构的指针，定义了文件的安全特性（如果操作系统支持的话）
								 OPEN_ALWAYS,	   /* CREATE_NEW: 创建文件, 如文件存在则会出错; CREATE_ALWAYS: 创建文件，会改写前一个文件; OPEN_EXISTING: 文件必须已经存在, 由设备提出要求;
												      OPEN_ALWAYS: 如文件不存在则创建它; TRUNCATE_EXISTING: 将现有文件缩短为零长度 */
						         flags,            // 文件属性
						         NULL);            // 如果不为零，则指定一个文件句柄。新文件将从这个文件中复制扩展属性

	if (mmHandle == INVALID_HANDLE_VALUE)  // 创建文件失败
	{
		error_code = GetLastError();
		cout << "创建mmf失败:" << error_code << endl;
	}
	else                                  // 创建mmf成功
	{
		DWORD high_size;
		DWORD file_size = GetFileSize(mmHandle, &high_size);    // 获取文件大小

		if (file_size == BAD_POS && (error_code = GetLastError()) != SUCCESS)   // 创建文件失败
		{
			CloseHandle(mmHandle);   //关闭内存映射文件
			cout << "error：" << error_code << endl;
		}

		//创建文件映射，如果要创建内存页面文件的映射，第一个参数设置为INVALID_HANDLE_VALUE
		mmfm = CreateFileMapping(mmHandle,                     /* 指向创建映射对象的文件句柄. 如果需要和物理文件关联, 要确保物理文件创建的时候的访问模式和由flProtect参数指定的"保护标识"匹配，
															   比如：物理文件只读， 内存映射需要读写就会发生错误。推荐将物理文件使用独占方式创建 */
										NULL,                  // 指向一个SECURITY_ATTRIBUTES结构体，决定是否被子进程继承。如果为NULL,句柄不能被继承，一般设置成NULL
										PAGE_READWRITE,        // 当文件映射的时候，需要此参数设置保护标识. PAGE_READWRITE: 同时有读和写的权限
										_mmfsizehigh,          // dwMaximumSizeHigh：文件映射对象的最大值的高位（DWORD）
										_mmfsizelow,           // dwMaximumSizeLow: 文件映射对象的最大值的低位（DWORD），如果dwMaximumSizeHigh为0，文件映射对象的最大值为由hFile标识的文件的大小
										_sharedname);          // 文件映射对象的名称

		error_code = GetLastError();
		if (SUCCESS != error_code)
		{
			cout << "createFileMapping error: " << error_code << endl;
		}
		else
		{
			if (mmfm == NULL)
			{
				if (mmHandle != INVALID_HANDLE_VALUE)
				{
					CloseHandle(mmHandle);
				}
			}
			else
			{
				DWORD view_access = FILE_MAP_ALL_ACCESS;      // 读和写的访问权限

				// 获得映射视图. 映射文件视图到调用进程的地址空间. 映射一个文件，让其指定的文件部分在调用进程的地址空间可见
				// 如果函数成功，返回值是映射视图的开始位置.如果函数失败，返回值为NULL，可以通过调用GetLastError函数获得详细的错误信息
				_mmfm_base_address = (float*)MapViewOfFile(mmfm,          // hFileMappingObject: 一个打开的映射文件对象的句柄，这个句柄可以由CreateFileMapping和OpenFileMapping函数返回
														 view_access,     // dwDesiredAccess: 指定访问文件视图的类型
														 0,               // dwFileOffsetHigh: 指定开始映射文件偏移量的高位
														 0,               // dwFileOffsetLow: 指定开始映射文件偏移量的低位
														 0);              // dwNumberOfBytesToMap: 指定需要映射的文件的字节数量，如果dwNumberOfBytesToMap为0，映射整个的文件

				if (_mmfm_base_address == NULL)
				{
					error_code = GetLastError();
					if (error_code != SUCCESS)
					{
						cout << "error code " << error_code << endl;
					}
				}
			}
		}
	}
}


void MMF::mmfClose(float * & _mmfm_base_address)
{
	//卸载映射
	UnmapViewOfFile(_mmfm_base_address);
	//关闭内存映射文件
	CloseHandle(mmfm);
	//关闭文件
	CloseHandle(mmHandle);
}