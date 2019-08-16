#ifndef CEVNS_HDF5_FUNCTIONS_H
#define CEVNS_HDF5_FUNCTIONS_H

#include <hdf5.h>
#include <iostream>

#define CONFLEN 300
#define STRLEN 20

  typedef struct{
     char param_key[CONFLEN];
     char param_value[CONFLEN];
   } run_info_t;

  typedef struct{
    int32_t event_id;
    unsigned int sensor_id;
    uint64_t time_bin;
    unsigned int charge;
  } sns_data_t;

  typedef struct{
    int32_t event_id;
    float x;
    float y;
    float z;
    float time;
    float energy;
    char label[STRLEN];
    int particle_id;
    int hit_id;
  } hit_info_t;

  typedef struct{
    int32_t event_id;
    int particle_id;
    char name[STRLEN];
    char primary;
    char mother_name[STRLEN];
    float initial_x;
    float initial_y;
    float initial_z;
    float initial_t;
    float final_x;
    float final_y;
    float final_z;
    float final_t;
    char initial_volume[STRLEN];
    char final_volume[STRLEN];
    float initial_momentum_x;
    float initial_momentum_y;
    float initial_momentum_z;
    float final_momentum_x;
    float final_momentum_y;
    float final_momentum_z;
    float kin_energy;
    char creator_proc[CONFLEN];
  } particle_info_t;

  typedef struct{
    unsigned int sensor_id;
    float x;
    float y;
    float z;
  } sns_pos_t;

  hsize_t createCevnsRunType();
  hsize_t createCevnsSensorDataType();
  hsize_t createCevnsHitInfoType();
  hsize_t createCevnsParticleInfoType();
  hsize_t createCevnsSensorPosType();

  hid_t createCevnsTable(hid_t group, std::string& table_name, hsize_t memtype);
  hid_t createCevnsGroup(hid_t file, std::string& groupName);

  void writeCevnsRun(run_info_t* runData, hid_t dataset, hid_t memtype, hsize_t counter);
  void writeCevnsSnsData(sns_data_t* snsData, hid_t dataset, hid_t memtype, hsize_t counter);
  void writeCevnsHit(hit_info_t* hitInfo, hid_t dataset, hid_t memtype, hsize_t counter);
  void writeCevnsParticle(particle_info_t* particleInfo, hid_t dataset, hid_t memtype, hsize_t counter);
  void writeCevnsSnsPos(sns_pos_t* snsPos, hid_t dataset, hid_t memtype, hsize_t counter);


#endif
