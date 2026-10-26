/*   ColDet - C++ 3D Collision Detection Library
 *   Copyright (C) 2000-2013   Amir Geva
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 * 
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA  02111-1307, USA.
 *
 * Any comments, questions and bug reports send to:
 *   amirgeva@gmail.com
 *
 * Or visit the home page: http://sourceforge.net/projects/coldet/
 */
#pragma once

#include "multiobject.h"
#include "cdmath3d.h"
#include "sweep_prune.h"
#include <vector>
#include <list>

__CD__BEGIN

typedef std::vector<CollisionDetails> details_vec;

class MultiObjectImpl : public MultiObjectSystem
{
protected:
  struct Model
  {
    Model() : col_model(0), updater(0), valid(false) {}
    bool              valid;
    CollisionModel3D* col_model;
    Vector3D          position;
    float             radius;
    const TransformUpdater* updater;
  };
  typedef std::vector<Model> model_vec;
  typedef std::list<int> int_list;

  model_vec   m_Models;
  int_list    m_FreeSlots;
  details_vec m_Results;

  /** Where a model is: from its updater, or at its position */
  static Matrix3D placement(const Model& m)
  {
    if (m.updater)
      return *reinterpret_cast<const Matrix3D*>(m.updater->update());
    Matrix3D t=Matrix3D::Identity;
    t.m[3][0]=m.position.x;
    t.m[3][1]=m.position.y;
    t.m[3][2]=m.position.z;
    return t;
  }

  static void place(Model& m)
  {
    m.col_model->setTransform(placement(m).m[0]);
  }

  bool exact_check(int id1, int id2)
  {
    bool rc=false;
    Model& m1=m_Models[id1];
    Model& m2=m_Models[id2];
    CollisionDetails& cd=m_Results.back();
    if (m1.col_model && m2.col_model)
    {
      // The second placement is passed to the test rather than set on the model: objects
      // may share a collision model, which can only hold one transform
      place(m1);
      Matrix3D t2=placement(m2);
      if (m1.col_model->collision(m2.col_model,-1,0,t2.m[0]))
      {
        m1.col_model->getCollidingTriangles(cd.t1,cd.t2,false);
        m1.col_model->getCollisionPoint(cd.point,false);
        rc=true;
      }
    }
    else
    if (m1.col_model || m2.col_model)
    {
      Model* m=&m2;
      CollisionModel3D* cm=m1.col_model;
      if (!cm) 
      { 
        cm=m2.col_model; 
        m=&m1; 
        place(m2);
      }
      else
        place(m1);
      if (cm->sphereCollision(m->position.get(),m->radius))
      {
        cm->getCollidingTriangles(cd.t1,cd.t2,false);
        cm->getCollisionPoint(cd.point,false);
        rc=true;
      }
    }
    else
    {
      // two sphere objects
      rc=SphereSphereCollision(m1.position.get(),m1.radius,m2.position.get(),m2.radius,cd.point);
    }
    return rc;
  }

  int allocate_model()
  {
    int id=-1;
    if (!m_FreeSlots.empty())
    {
      id=m_FreeSlots.front();
      m_FreeSlots.pop_front();
    }
    else
    {
      id=int(m_Models.size());
      m_Models.push_back(Model());
    }
    return id;
  }

  /** Adds a model to the system's structure, or releases its slot if there is no room.
      Returns the id, or -1 */
  int insert(int id, const float* position, float radius)
  {
    if (add_object(id,position,radius)) return id;
    m_Models[id]=Model();
    if (id==int(m_Models.size())-1) m_Models.pop_back();
    else m_FreeSlots.push_front(id);
    return -1;
  }
public:
  MultiObjectImpl()
  {}

  virtual int add_object(CollisionModel3D* model, const TransformUpdater* updater)
  {
    int id=-1;
    if (!model) return id;
    id=allocate_model();
    Model& m=m_Models[id];
    m.col_model=model;
    m.valid=true;
    m.radius=model->getRadius();
    m.updater=updater;
    const float* mat=updater->update();
    m.position.x=mat[12];
    m.position.y=mat[13];
    m.position.z=mat[14];
    return insert(id,&mat[12],m.radius);
  }

  virtual int add_object(CollisionModel3D* model, const float* position)
  {
    int id=-1;
    if (!model) return id;
    id=allocate_model();
    Model& m=m_Models[id];
    m.col_model=model;
    m.valid=true;
    m.radius=model->getRadius();
    m.updater=0;
    m.position.x=position[0];
    m.position.y=position[1];
    m.position.z=position[2];
    return insert(id,position,m.radius);
  }

  virtual int add_sphere_object(float radius, const float* position)
  {
    int id=-1;
    id=allocate_model();
    Model& m=m_Models[id];
    m.valid=true;
    m.radius=radius;
    m.updater=0;
    m.position.x=position[0];
    m.position.y=position[1];
    m.position.z=position[2];
    return insert(id,position,radius);
  }

  /** Returns false if the system has no room for the object */
  virtual bool add_object(int id, const float* position, float radius) = 0;
  virtual void remove_object_from_system(int id) = 0;

  virtual void remove_object(int id)
  {
    if (id<0 || id>=int(m_Models.size())) return;
    Model& m=m_Models[id];
    if (!m.valid) return;
    remove_object_from_system(id);
    m.valid=false;
    m.col_model=0;
    m.position=Vector3D(0,0,0);
    m.updater=0;
    m_FreeSlots.push_back(id);
  }

  CollisionModel3D* get_collision_model(int id)
  {
    if (id<0 || id>=int(m_Models.size())) return 0;
    Model& m=m_Models[id];
    if (!m.valid) return 0;
    return m.col_model;
  }

  virtual void move_object(int id, const float* new_position)
  {
    if (id<0 || id>=int(m_Models.size())) return;
    Model& m=m_Models[id];
    if (!m.valid) return;
    Vector3D pos(new_position[0],new_position[1],new_position[2]);
    Vector3D offset=pos-m.position;
    update_position(id,m,pos,offset);
    m.position=pos;
  }

  virtual int check_collisions(bool exact)
  {
    m_Results.clear();
    return 0;
  }

  virtual bool get_collision(int index, CollisionDetails& details)
  {
    if (index<0 || index>=int(m_Results.size())) return false;
    CollisionDetails& cd=m_Results[index];
    details = cd;
    return true;
  }

  virtual void update_position(int id, Model& m, const Vector3D& pos, const Vector3D& delta)
  {}

  /** Objects placed by an updater move with it: take up their current positions */
  void follow_updaters()
  {
    for(int id=0;id<int(m_Models.size());++id)
    {
      Model& m=m_Models[id];
      if (!m.valid || !m.updater) continue;
      const float* mat=m.updater->update();
      Vector3D pos(mat[12],mat[13],mat[14]);
      if (pos.x!=m.position.x || pos.y!=m.position.y || pos.z!=m.position.z)
      {
        update_position(id,m,pos,pos-m.position);
        m.position=pos;
      }
    }
  }
};

class SweepPruneSystem : public MultiObjectImpl
{
  SNP::SweepPrune m_SP;
public:
  SweepPruneSystem(int max_objects) : m_SP(max_objects) {}

  virtual void remove_object_from_system(int id)
  {
    m_SP.remove_box(id);
  }

  virtual bool add_object(int id, const float* position, float radius)
  {
    return m_SP.add_box(id,position,radius);
  }

  virtual void update_position(int id, Model& m, const Vector3D& pos, const Vector3D& delta)
  {
    m_SP.move_box(id,delta.get());
  }

  virtual int check_collisions(bool exact)
  {
    follow_updaters();
    m_Results.clear();
    const SNP::intr_set& s=m_SP.get_intersections();
    for(SNP::intr_set::const_iterator it=s.begin();it!=s.end();++it)
    {
	  const SNP::Intersection& intr=*it;
      m_Results.push_back(CollisionDetails());
      CollisionDetails& cd=m_Results.back();
      cd.id1=intr.id1;
      cd.id2=intr.id2;
      if (exact)
      {
        if (!exact_check(intr.id1,intr.id2)) m_Results.pop_back();
      }
    }
    return m_Results.size();
  }
};

class SphereSystem : public MultiObjectImpl
{
public:
  virtual bool add_object(int id, const float* position, float radius) { return true; }
  virtual void remove_object_from_system(int id) {}
  virtual int check_collisions(bool exact)
  {
    follow_updaters();
    m_Results.clear();
    int n=m_Models.size();
    for(int id1=0;id1<n;++id1)
    {
      Model& m1=m_Models[id1];
      if (!m1.valid) continue;
      float r1=m1.radius;
      for(int id2=id1+1;id2<n;++id2)
      {
        Model& m2=m_Models[id2];
        if (!m2.valid) continue;
        float r2=m2.radius;
        float sumr=r1+r2;
        float sqd=(m1.position-m2.position).SquareMagnitude();
        if ((sumr*sumr)>sqd)
        {
          m_Results.push_back(CollisionDetails());
          CollisionDetails& cd=m_Results.back();
          cd.id1=id1;
          cd.id2=id2;
          if (exact) 
            if (!exact_check(id1,id2)) m_Results.pop_back();
        }
      }
    }
    return m_Results.size();
  }


};

__CD__END


