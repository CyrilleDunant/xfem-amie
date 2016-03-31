//
// C++ Implementation: physics_base
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "physics_base.h"

using namespace Amie ;


LinearForm::LinearForm(const Matrix & p, bool t, bool s, size_t numdof ) : Form(p, t, s, numdof)
{
    this->type = PURE_LINEAR ;
}

LinearForm::LinearForm(bool t, bool s, size_t numdof ) : Form(Matrix(), t, s, numdof)
{
    this->type = PURE_LINEAR ;
}

LinearForm::~LinearForm() { }

bool LinearForm::fractured() const {
    return false  ;
} 

void LinearForm::step(double timestep, ElementState & s, double maxscore)
{

}


NonLinearForm::NonLinearForm() : Form(Matrix(), false, false, 2)
{
    this->type = NON_LINEAR ;
}

NonLinearForm::~NonLinearForm() { } 

std::vector<Point> NonLinearForm::getIntegrationHints()
{
    return hints ;
}

Point NonLinearForm::getIntegrationHint(size_t i)
{
    return hints[i] ;
}

void NonLinearForm::setIntegrationHints(std::vector<Point> h)
{
    hints = h ;
}

void NonLinearForm::step(double timestep, ElementState & s, double maxscore)
{
}

bool NonLinearForm::fractured() const
{
    return false ;
}


LinearFormAndConstant::LinearFormAndConstant(const Matrix & p, Matrix c) : Form(p, false, false, p.numRows()/3+1)
{
    this->type = LINEAR_AND_CONSTANT ;
}

LinearFormAndConstant::~LinearFormAndConstant() { } 

