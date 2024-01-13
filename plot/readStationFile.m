function [locations, stationNames] = readStationFile(file_path)
    % 读取地震台站文本文件
    fid = fopen(file_path, 'r');
    
    % 检查文件是否成功打开
    if fid == -1
        error('无法打开文件。');
    end
    
    % 初始化存储数组
    data = [];
    stationNames = {};
    
    % 逐行读取文件内容
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)
            % 使用空格或逗号拆分行
            parts = strsplit(line, {' ', ','});
            
            % 提取经度、纬度、深度信息
            longitude = str2double(parts{1});
            latitude = str2double(parts{2});
            depth = str2double(parts{3});
            
            % 提取台网名和台站名
            networkName = parts{4};
            stationName = parts{5};
            
            % 存储位置信息
            data = [data; longitude, latitude, depth];
            
            % 存储台站名，使用'.'连接台网名和台站名
            stationNames = [stationNames; strcat(networkName, '.', stationName)];
        end
    end
    
    % 关闭文件
    fclose(fid);
    
    % 返回位置信息和台站名
    locations = data;
end
