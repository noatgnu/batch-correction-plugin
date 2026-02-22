# Batch Effect Correction

**ID**: `batch-correction`  
**Version**: 1.0.0  
**Category**: preprocessing  
**Author**: CauldronGO Team

## Description

Correct for batch effects in proteomics data

## Runtime

- **Type**: `r`
- **Script**: `batch_correction.R`

## Inputs

| Name | Label | Type | Required | Default | Visibility |
|------|-------|------|----------|---------|------------|
| `input_file` | Input File | file | Yes | - | Always visible |
| `annotation_file` | Annotation File | file | Yes | - | Always visible |
| `columns_name` | Sample Columns | column-selector (multiple) | Yes | - | Always visible |
| `method` | Correction Method | select (combat, limma, ruvseq) | Yes | combat | Always visible |
| `preserve_column` | Preserve Column | text | No | - | Always visible |
| `use_log2` | Use Log2 | boolean | No | false | Always visible |

### Input Details

#### Input File (`input_file`)

Data file with batch effects


#### Annotation File (`annotation_file`)

Sample annotation file with 'Sample' and 'Batch' columns


#### Sample Columns (`columns_name`)

Select columns containing sample data

- **Column Source**: `input_file`

#### Correction Method (`method`)

Method to use for batch correction

- **Options**: `combat`, `limma`, `ruvseq`

#### Preserve Column (`preserve_column`)

Column name from annotation file to preserve during correction (e.g., Condition, BioReplicate)


#### Use Log2 (`use_log2`)

Apply log2 transformation before correction


## Outputs

| Name | File | Type | Format | Description |
|------|------|------|--------|-------------|
| `corrected_data` | `batch_corrected.data.txt` | data | tsv | Batch-corrected data matrix |
| `batch_info` | `batch_info.txt` | data | tsv | Batch correction summary information |

## Sample Annotation

This plugin supports sample annotation:

- **Samples From**: `columns_name`
- **Annotation File**: `annotation_file`

## Requirements

- **R**: >=4.0
- **Packages**:
  - sva
  - limma
  - RUVSeq

## Example Data

This plugin includes example data for testing:

```yaml
  columns_name: [C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_BCA_LT-IP_01.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_BCA_LT-IP_02.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_BCA_LT-MockIP_01.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_BCA_LT-MockIP_02.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_Pepide-CBQCA_LT-IP_01.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_Pepide-CBQCA_LT-IP_02.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_Pepide-CBQCA_LT-MockIP_01.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_Pepide-CBQCA_LT-MockIP_02.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_BCA_LT-WCL_01.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_BCA_LT-WCL_02.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_Pepide-CBQCA_LT-WCL_01.raw C:\Raja\DIA-NN searches\June 2022\LT-CBQCA-Test_DIA\RN-DS_220106_Pepide-CBQCA_LT-WCL_02.raw]
  method: combat
  use_log2: false
  input_file: diann/imputed.data.txt
  annotation_file: differential_analysis/batch_info.txt
  columns_name_source: diann/imputed.data.txt
```

Load example data by clicking the **Load Example** button in the UI.

## Usage

### Via UI

1. Navigate to **preprocessing** → **Batch Effect Correction**
2. Fill in the required inputs
3. Click **Run Analysis**

### Via Plugin System

```typescript
const jobId = await pluginService.executePlugin('batch-correction', {
  // Add parameters here
});
```
