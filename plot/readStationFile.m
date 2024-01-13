function [locations, stationNames] = readStationFile(file_path)
    % ��ȡ����̨վ�ı��ļ�
    fid = fopen(file_path, 'r');
    
    % ����ļ��Ƿ�ɹ���
    if fid == -1
        error('�޷����ļ���');
    end
    
    % ��ʼ���洢����
    data = [];
    stationNames = {};
    
    % ���ж�ȡ�ļ�����
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)
            % ʹ�ÿո�򶺺Ų����
            parts = strsplit(line, {' ', ','});
            
            % ��ȡ���ȡ�γ�ȡ������Ϣ
            longitude = str2double(parts{1});
            latitude = str2double(parts{2});
            depth = str2double(parts{3});
            
            % ��ȡ̨������̨վ��
            networkName = parts{4};
            stationName = parts{5};
            
            % �洢λ����Ϣ
            data = [data; longitude, latitude, depth];
            
            % �洢̨վ����ʹ��'.'����̨������̨վ��
            stationNames = [stationNames; strcat(networkName, '.', stationName)];
        end
    end
    
    % �ر��ļ�
    fclose(fid);
    
    % ����λ����Ϣ��̨վ��
    locations = data;
end
