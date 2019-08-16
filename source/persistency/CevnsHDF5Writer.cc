#include "CevnsHDF5Writer.h"

#include <sstream>
#include <cstring>
#include <stdlib.h>
#include <vector>

#include <stdint.h>
#include <iostream>

using namespace nexus;


CevnsHDF5Writer::CevnsHDF5Writer():
  _file(0), _irun(0), _ismp(0), _ihit(0),
  _ipart(0), _ipos(0)
{
}

CevnsHDF5Writer::~CevnsHDF5Writer()
{
}

void CevnsHDF5Writer::Open(std::string fileName)
{
  _firstEvent= true;

  _file = H5Fcreate( fileName.c_str(), H5F_ACC_TRUNC,
                      H5P_DEFAULT, H5P_DEFAULT );

  std::string group_name = "/MC";
  _group = createCevnsGroup(_file, group_name);

  std::string run_table_name = "configuration";
  _memtypeRun = createCevnsRunType();
  _runTable = createCevnsTable(_group, run_table_name, _memtypeRun);

  std::string sns_data_table_name = "waveforms";
  _memtypeSnsData = createCevnsSensorDataType();
  _snsDataTable = createCevnsTable(_group, sns_data_table_name, _memtypeSnsData);

  std::string hit_info_table_name = "hits";
  _memtypeHitInfo = createCevnsHitInfoType();
  _hitInfoTable = createCevnsTable(_group, hit_info_table_name, _memtypeHitInfo);

  std::string particle_info_table_name = "particles";
  _memtypeParticleInfo = createCevnsParticleInfoType();
  _particleInfoTable = createCevnsTable(_group, particle_info_table_name, _memtypeParticleInfo);

  std::string sns_pos_table_name = "sensor_positions";
  _memtypeSnsPos = createCevnsSensorPosType();
  _snsPosTable = createCevnsTable(_group, sns_pos_table_name, _memtypeSnsPos);

  _isOpen = true;
}

void CevnsHDF5Writer::Close()
{
  _isOpen=false;
  H5Fclose(_file);
}

void CevnsHDF5Writer::WriteRunInfo(const char* param_key, const char* param_value)
{
  run_info_t runData;
  memset(runData.param_key,   0, CONFLEN);
  memset(runData.param_value, 0, CONFLEN);
  strcpy(runData.param_key, param_key);
  strcpy(runData.param_value, param_value);
  writeCevnsRun(&runData, _runTable, _memtypeRun, _irun);

  _irun++;
}

void CevnsHDF5Writer::WriteSensorDataInfo(int evt_number, unsigned int sensor_id, unsigned int time_bin, unsigned int charge)
{
  sns_data_t snsData;
  snsData.event_id = evt_number;
  snsData.sensor_id = sensor_id;
  snsData.time_bin = time_bin;
  snsData.charge = charge;
  writeCevnsSnsData(&snsData, _snsDataTable, _memtypeSnsData, _ismp);

  _ismp++;
}


void CevnsHDF5Writer::WriteHitInfo(int evt_number, int particle_indx, int hit_indx, float hit_position_x, float hit_position_y, float hit_position_z, float hit_time, float hit_energy, const char* label)
{
  hit_info_t trueInfo;
  trueInfo.event_id = evt_number;
  memset(trueInfo.label, 0, STRLEN);
  //memcpy(trueInfo.hit_position, hit_position, sizeof(*hit_position)*size_position);
  trueInfo.x = hit_position_x;
  trueInfo.y = hit_position_y;
  trueInfo.z = hit_position_z;
  trueInfo.time = hit_time;
  trueInfo.energy = hit_energy;
  strcpy(trueInfo.label, label);
  trueInfo.particle_id = particle_indx;
  trueInfo.hit_id = hit_indx;
  writeCevnsHit(&trueInfo,  _hitInfoTable, _memtypeHitInfo, _ihit);

  _ihit++;
}

void CevnsHDF5Writer::WriteParticleInfo(int evt_number, int particle_indx, const char* particle_name, char primary, const char* mother_name, float initial_vertex_x, float initial_vertex_y, float initial_vertex_z, float initial_vertex_t, float final_vertex_x, float final_vertex_y, float final_vertex_z, float final_vertex_t, const char* initial_volume, const char* final_volume, float momentum_x, float momentum_y, float momentum_z, float final_momentum_x, float final_momentum_y, float final_momentum_z, float kin_energy, const char* creator_proc)
{
  particle_info_t trueInfo;
  trueInfo.event_id = evt_number;
  trueInfo.particle_id = particle_indx;
  memset(trueInfo.name, 0, STRLEN);
  strcpy(trueInfo.name, particle_name);
  trueInfo.primary = primary;
  memset(trueInfo.mother_name, 0, STRLEN);
  strcpy(trueInfo.mother_name, mother_name);
  trueInfo.initial_x = initial_vertex_x;
  trueInfo.initial_y = initial_vertex_y;
  trueInfo.initial_z = initial_vertex_z;
  trueInfo.initial_t = initial_vertex_t;
  trueInfo.final_x = final_vertex_x;
  trueInfo.final_y = final_vertex_y;
  trueInfo.final_z = final_vertex_z;
  trueInfo.final_t = final_vertex_t;
  // memcpy(trueInfo.initial_vertex, initial_vertex, sizeof(*initial_vertex)*size_initial_vertex);
  //  memcpy(trueInfo.final_vertex, final_vertex, sizeof(*final_vertex)*size_final_vertex);
  memset(trueInfo.initial_volume, 0, STRLEN);
  strcpy(trueInfo.initial_volume, initial_volume);
  memset(trueInfo.final_volume, 0, STRLEN);
  strcpy(trueInfo.final_volume, final_volume);
  trueInfo.initial_momentum_x = momentum_x;
  trueInfo.initial_momentum_y = momentum_y;
  trueInfo.initial_momentum_z = momentum_z;
  trueInfo.final_momentum_x = final_momentum_x;
  trueInfo.final_momentum_y = final_momentum_y;
  trueInfo.final_momentum_z = final_momentum_z;
  //  memcpy(trueInfo.momentum, momentum, sizeof(*momentum)*size_momentum);
  trueInfo.kin_energy = kin_energy;
  memset(trueInfo.creator_proc, 0, STRLEN);
  strcpy(trueInfo.creator_proc, creator_proc);
  writeCevnsParticle(&trueInfo,  _particleInfoTable, _memtypeParticleInfo, _ipart);

  _ipart++;
}

void CevnsHDF5Writer::WriteSensorPosInfo(unsigned int sensor_id, float x, float y, float z)
{
  sns_pos_t snsPos;
  snsPos.sensor_id = sensor_id;
  snsPos.x = x;
  snsPos.y = y;
  snsPos.z = z;
  writeCevnsSnsPos(&snsPos, _snsPosTable, _memtypeSnsPos, _ipos);

  _ipos++;
}
