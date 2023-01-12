
typedef struct
{
	double time;         //time (s)
	double wb_ib[3];     //gyro data(deg/s)
	double fb_ib[3];     //acclerometer(g)
	double temperature;   //temperature (C)
}IMUData_t;

typedef struct
{
	double temperature;
	double Ka[9];       //acc calibration matrix (列优先存储,类似matlab)
	double Kg[9];       //gyro calibration matrix(列优先存储,类似matlab)
	double ba[3];       //bias acc	(g)
	double bg[3];       //bias gyro	(deg/s)
}IMUParam_t;

typedef struct
{
#define kMAXTEMPNUM 10
	int size;
	IMUParam_t imu_param[kMAXTEMPNUM];
}AllTempIMUParam_t;

// 线性插值
bool ParamInterp(const AllTempIMUParam_t &all_param, const double &temperature, IMUParam_t &interp_param)
{
	if (all_param.size < 1) return false;

	if (all_param.size == 1)
	{
		interp_param = all_param.imu_param[0];
		return true;
	}

	double temp = temperature;
	double ba_x = 0.0, ba_y = 0.0, ba_z = 0.0;
	double bg_x = 0.0, bg_y = 0.0, bg_z = 0.0;

	double Kaxx = 0.0, Kaxy = 0.0, Kaxz = 0.0;
	double Kayx = 0.0, Kayy = 0.0, Kayz = 0.0;
	double Kazx = 0.0, Kazy = 0.0, Kazz = 0.0;

	double Kgxx = 0.0, Kgxy = 0.0, Kgxz = 0.0;	//036
	double Kgyx = 0.0, Kgyy = 0.0, Kgyz = 0.0;	//147
	double Kgzx = 0.0, Kgzy = 0.0, Kgzz = 0.0;	//258

	//从全温参数中找到最接近检校温度的两个温度点
	double max_temp = -99999.9, min_temp = 99999.9;
	int temp_1 = 0, temp_2 = 0;

	for (int i = 0; i < all_param.size; i++)

	{
		if (all_param.imu_param[i].temperature > max_temp)
		{
			max_temp = all_param.imu_param[i].temperature;
			temp_2 = i;
		}
		if (all_param.imu_param[i].temperature < min_temp)
		{
			min_temp = all_param.imu_param[i].temperature;
			temp_1 = i;
		}
	}

	//如果在区间内
	if (temp >= min_temp&&temp <= max_temp)
	{
		double upper = 999.0, lower = -999.0;
		for (int i = 0; i < all_param.size; i++)
		{
			if (all_param.imu_param[i].temperature>temp&& all_param.imu_param[i].temperature<all_param.imu_param[temp_2].temperature)
			{
				temp_2 = i;
			}
			if (all_param.imu_param[i].temperature<temp&& all_param.imu_param[i].temperature>all_param.imu_param[temp_1].temperature)
			{
				temp_1 = i;
			}
		}
	}

	//如果在区间外
	if (temp<min_temp)
	{
		temp_2 = (temp_1 + 1) % all_param.size;
		for (int i = 0; i < all_param.size; i++)
		{
			if (all_param.imu_param[i].temperature<all_param.imu_param[temp_2].temperature && i != temp_1)
			{
				temp_2 = i;
			}
		}
	}

	//如果在区间外
	if (temp>max_temp)
	{
		temp_1 = (temp_2 + 1) % all_param.size;
		for (int i = 0; i < all_param.size; i++)
		{
			if (all_param.imu_param[i].temperature>all_param.imu_param[temp_1].temperature&& i != temp_2)
			{
				temp_1 = i;
			}
		}
	}

	///开始插值		// y=k1*y1+k2*y2
	double k1;		// k1=(x-x2)/(x1-x2)
	double k2;		// k2=(x-x1)/(x2-x1)
	k1 = (temp - all_param.imu_param[temp_2].temperature) / (all_param.imu_param[temp_1].temperature - all_param.imu_param[temp_2].temperature);
	k2 = (temp - all_param.imu_param[temp_1].temperature) / (all_param.imu_param[temp_2].temperature - all_param.imu_param[temp_1].temperature);

	ba_x = k1*all_param.imu_param[temp_1].ba[0] + k2*all_param.imu_param[temp_2].ba[0];
	ba_y = k1*all_param.imu_param[temp_1].ba[1] + k2*all_param.imu_param[temp_2].ba[1];
	ba_z = k1*all_param.imu_param[temp_1].ba[2] + k2*all_param.imu_param[temp_2].ba[2];

	bg_x = k1*all_param.imu_param[temp_1].bg[0] + k2*all_param.imu_param[temp_2].bg[0];
	bg_y = k1*all_param.imu_param[temp_1].bg[1] + k2*all_param.imu_param[temp_2].bg[1];
	bg_z = k1*all_param.imu_param[temp_1].bg[2] + k2*all_param.imu_param[temp_2].bg[2];

	Kaxx = k1*all_param.imu_param[temp_1].Ka[0] + k2*all_param.imu_param[temp_2].Ka[0];
	Kayx = k1*all_param.imu_param[temp_1].Ka[1] + k2*all_param.imu_param[temp_2].Ka[1];
	Kazx = k1*all_param.imu_param[temp_1].Ka[2] + k2*all_param.imu_param[temp_2].Ka[2];
	Kaxy = k1*all_param.imu_param[temp_1].Ka[3] + k2*all_param.imu_param[temp_2].Ka[3];
	Kayy = k1*all_param.imu_param[temp_1].Ka[4] + k2*all_param.imu_param[temp_2].Ka[4];
	Kazy = k1*all_param.imu_param[temp_1].Ka[5] + k2*all_param.imu_param[temp_2].Ka[5];
	Kaxz = k1*all_param.imu_param[temp_1].Ka[6] + k2*all_param.imu_param[temp_2].Ka[6];
	Kayz = k1*all_param.imu_param[temp_1].Ka[7] + k2*all_param.imu_param[temp_2].Ka[7];
	Kazz = k1*all_param.imu_param[temp_1].Ka[8] + k2*all_param.imu_param[temp_2].Ka[8];

	Kgxx = k1*all_param.imu_param[temp_1].Kg[0] + k2*all_param.imu_param[temp_2].Kg[0];
	Kgyx = k1*all_param.imu_param[temp_1].Kg[1] + k2*all_param.imu_param[temp_2].Kg[1];
	Kgzx = k1*all_param.imu_param[temp_1].Kg[2] + k2*all_param.imu_param[temp_2].Kg[2];
	Kgxy = k1*all_param.imu_param[temp_1].Kg[3] + k2*all_param.imu_param[temp_2].Kg[3];
	Kgyy = k1*all_param.imu_param[temp_1].Kg[4] + k2*all_param.imu_param[temp_2].Kg[4];
	Kgzy = k1*all_param.imu_param[temp_1].Kg[5] + k2*all_param.imu_param[temp_2].Kg[5];
	Kgxz = k1*all_param.imu_param[temp_1].Kg[6] + k2*all_param.imu_param[temp_2].Kg[6];
	Kgyz = k1*all_param.imu_param[temp_1].Kg[7] + k2*all_param.imu_param[temp_2].Kg[7];
	Kgzz = k1*all_param.imu_param[temp_1].Kg[8] + k2*all_param.imu_param[temp_2].Kg[8];
	interp_param.temperature = temperature;

	interp_param.ba[0] = ba_x;
	interp_param.ba[1] = ba_y;
	interp_param.ba[2] = ba_z;

	interp_param.bg[0] = bg_x;
	interp_param.bg[1] = bg_y;
	interp_param.bg[2] = bg_z;

	interp_param.Ka[0] = Kaxx;	interp_param.Ka[3] = Kaxy;	interp_param.Ka[6] = Kaxz;
	interp_param.Ka[1] = Kayx;	interp_param.Ka[4] = Kayy;	interp_param.Ka[7] = Kayz;
	interp_param.Ka[2] = Kazx;	interp_param.Ka[5] = Kazy;	interp_param.Ka[8] = Kazz;

	interp_param.Kg[0] = Kgxx;	interp_param.Kg[3] = Kgxy;	interp_param.Kg[6] = Kgxz;
	interp_param.Kg[1] = Kgyx;	interp_param.Kg[4] = Kgyy;	interp_param.Kg[7] = Kgyz;
	interp_param.Kg[2] = Kgzx;	interp_param.Kg[5] = Kgzy;	interp_param.Kg[8] = Kgzz;

	return true;
}

// 获取标定参数补偿后的IMU数据
void UpdateCompensatedIMUData(const AllTempIMUParam_t &all_param, IMUData_t &imu_data)
{	
	IMUParam_t interp_param = { 0 };
	IMUData_t dateTemp = imu_data;

	if (ParamInterp(all_param, imu_data.temperature, interp_param))
	{
		dateTemp.fb_ib[0] -= interp_param.ba[0]; //请确认零偏的单位与imu数据的单位一致。
		dateTemp.fb_ib[1] -= interp_param.ba[1];
		dateTemp.fb_ib[2] -= interp_param.ba[2];

		dateTemp.wb_ib[0] -= interp_param.bg[0];
		dateTemp.wb_ib[1] -= interp_param.bg[1];
		dateTemp.wb_ib[2] -= interp_param.bg[2];

		imu_data.fb_ib[0] = interp_param.Ka[0] * dateTemp.fb_ib[0] + interp_param.Ka[3] * dateTemp.fb_ib[1] + interp_param.Ka[6] * dateTemp.fb_ib[2];
		imu_data.fb_ib[1] = interp_param.Ka[1] * dateTemp.fb_ib[0] + interp_param.Ka[4] * dateTemp.fb_ib[1] + interp_param.Ka[7] * dateTemp.fb_ib[2];
		imu_data.fb_ib[2] = interp_param.Ka[2] * dateTemp.fb_ib[0] + interp_param.Ka[5] * dateTemp.fb_ib[1] + interp_param.Ka[8] * dateTemp.fb_ib[2];

		imu_data.wb_ib[0] = interp_param.Kg[0] * dateTemp.wb_ib[0] + interp_param.Kg[3] * dateTemp.wb_ib[1] + interp_param.Kg[6] * dateTemp.wb_ib[2];
		imu_data.wb_ib[1] = interp_param.Kg[1] * dateTemp.wb_ib[0] + interp_param.Kg[4] * dateTemp.wb_ib[1] + interp_param.Kg[7] * dateTemp.wb_ib[2];
		imu_data.wb_ib[2] = interp_param.Kg[2] * dateTemp.wb_ib[0] + interp_param.Kg[5] * dateTemp.wb_ib[1] + interp_param.Kg[8] * dateTemp.wb_ib[2];
	}
}



