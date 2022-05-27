USE [W101_Covid19_EpiCell_Workspace]
GO


-- INDEX PUSHED TABLE
---- first delete existing index
DROP INDEX IF EXISTS [nhs_index] on [dbo].[Severity_BA1_BA2]; -- clustered index on NHS number
DROP INDEX IF EXISTS [finalid_index] on [dbo].[Severity_BA1_BA2]; -- non-clustered index on Final ID

---- now create new index 
CREATE CLUSTERED INDEX [nhs_index] on [dbo].[Severity_BA1_BA2]
([nhs_number] ASC
) WITH (STATISTICS_NORECOMPUTE = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY];  

CREATE NONCLUSTERED INDEX [finalid_index] on [dbo].[Severity_BA1_BA2]
([final_id] ASC
) WITH (STATISTICS_NORECOMPUTE = OFF, SORT_IN_TEMPDB = OFF, IGNORE_DUP_KEY = OFF, ONLINE = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY];  


-- LINKAGE
IF OBJECT_ID('[dbo].[Severity_BA1_BA2_NIMS]','U') IS NOT NULL 
DROP TABLE [W101_Covid19_EpiCell_Workspace].[dbo].[Severity_BA1_BA2_NIMS]; -- CHECK IF TABLE EXISTS, IF SO THEN DELETE

SELECT D.*
	  ,NIMS.*
	  ,dbo.fnNhsNumberClassification(D.nhs_number) AS NHS_Check
		,CASE WHEN NIMS.NHSNumber IS NULL THEN 'UNLINKED'
		WHEN DATE_D1 IS NOT NULL OR DATE_D2 IS NOT NULL OR DATE_D3 IS NOT NULL THEN 'VACCINATED'
		ELSE 'UNVACCINATED'
	    END AS v_status

INTO [W101_Covid19_EpiCell_Workspace].[dbo].[Severity_BA1_BA2_NIMS]

FROM [W101_Covid19_EpiCell_Workspace].[dbo].[Severity_BA1_BA2] D 

LEFT JOIN 
( SELECT

	   N.NHSNumber
			,MAX(CASE WHEN N.dv_dose='Dose1' THEN N.dv_dose END) AS DOSE1
			,MAX(CASE WHEN N.dv_dose='Dose2' THEN N.dv_dose END) AS DOSE2
			,MAX(CASE WHEN N.dv_dose='Dose3' OR N.dv_dose='Booster1' THEN N.dv_dose END) AS DOSE3
			,MAX(CASE WHEN N.dv_dose='Dose1' THEN N.dv_DateAdministered END) AS DATE_D1
			,MAX(CASE WHEN N.dv_dose='Dose2' THEN N.dv_DateAdministered END) AS DATE_D2
			,MAX(CASE WHEN N.dv_dose='Dose3' OR N.dv_dose='Booster1' THEN N.dv_DateAdministered END) AS DATE_D3


FROM [W101_Covid19_EpiCell_Workspace].[dbo].[v_W101_VaccineData_Default_COVID19_Vaccinations] N

GROUP BY N.NHSNumber) AS NIMS

ON D.nhs_number = NIMS.NHSNumber