/** @file differencefinite.h */

#ifndef DIFFFINITE_H
#define DIFFFINITE_H

#include <vector>
#include<algorithm>

#include "toolbox.h"

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
  * @brief Énumérateur sur les schémas DF disponibles
  */
typedef enum
{
    ORDER_1_FORWARD,
    ORDER_2_FORWARD,
    ORDER_3_FORWARD,
    ORDER_4_FORWARD,
    ORDER_5_FORWARD,
    ORDER_6_FORWARD,
    ORDER_7_FORWARD,
    ORDER_8_FORWARD,
    ORDER_1_CENTRAL,
    ORDER_2_CENTRAL,
    ORDER_3_CENTRAL,
    ORDER_4_CENTRAL,
    ORDER_5_CENTRAL,
    ORDER_6_CENTRAL,
    ORDER_7_CENTRAL,
    ORDER_8_CENTRAL,
    ORDER_1_BACKWARD,
    ORDER_2_BACKWARD,
    ORDER_3_BACKWARD,
    ORDER_4_BACKWARD,
    ORDER_5_BACKWARD,
    ORDER_6_BACKWARD,
    ORDER_7_BACKWARD,
    ORDER_8_BACKWARD,
} ORDERS;


/**
 * @brief Structure DF regroupant les indices des points du schémas et les poids associés
 */
struct DFIdxCoeff
{
    std::vector<int> idxs = {};
    std::vector<double> coeffs = {};

};

/**
 * @brief Structure regroupant les ordres disponibles (voir DFIdxCoeff).
 */
struct DFOrders
{
    DFIdxCoeff Order1;
    DFIdxCoeff Order2;
    DFIdxCoeff Order3;
    DFIdxCoeff Order4;
    DFIdxCoeff Order5;
    DFIdxCoeff Order6;
    DFIdxCoeff Order7;
    DFIdxCoeff Order8;
};

/**
 * @brief Structure DF particulière pour la dérivée n-ième.
 */
struct DerivativeStruct
{
    DFOrders Backward;
    DFOrders Central;
    DFOrders Forward;
};

/**
 * @brief Structure DF regroupant les dérivées actuellement disponibles.
 */
struct DFStruct
{
    /**
     * @brief Dérivée première
     */
    DerivativeStruct Derivative_1;

    /**
     * @brief Dérivée seconde
     */
    DerivativeStruct Derivative_2;

    /**
     * @brief Cosntructeur
     */
    DFStruct();
};

/**
 * @brief Constructeur de schéma DF.
 * @param der degré de la dérivée à approximer.
 * @param order ordre demandé avec précision amont, aval ou centré
 * @param idxs pointeur vers le stockage des indices de points.
 * @param coeffs pointeur vers le stockage des poids.
 */
void DFOrderBuild (int der, ORDERS order, std::vector<int>* idxs, std::vector<double>* coeffs);

/**
 * @brief Inversion d'une structure DF (voir le lien entre avant et aval).
 * @param O structure à inverser.
 * @return la structure O inversé.
 */
DFOrders Reverse (DFOrders O);

/**
 * @brief Inversion d'une structure DF (voir le lien entre avant et aval).
 * @param D structure à inverser.
 * @return la structure D inversé.
 */
DFIdxCoeff Reverse (DFIdxCoeff D);

/** @} */
#endif // DIFFFINITE_H
