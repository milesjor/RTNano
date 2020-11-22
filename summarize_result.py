import pandas as pd
import logging


def rule():
    rules = {
        'POS+++': '3 region_alignment_record >= 50',
        'POS++': '2 region_alignment_record >= 50',
        'POS+': '1 region_alignment_record >= 50'
                'OR 2 region_alignment_record >= 20'
                'OR 3 region_alignment_record >= 10'
                'OR 4 region_alignment_record >= 5',
        'POS-': 'ACTB house_keeping_record >= 1000'
                'AND any region_alignment_record >= 1',
        'NEG': 'ACTB house_keeping_record >= 1000'
                'AND all region_alignment_record = 0',
        'UNK': 'ACTB house_keeping_record < 1000'
                'AND all region_alignment_record = 0',
    }
    return rules


def summarize(updated_result, target_amplicon, housekeeping, NTC):
    # updated_result = '/Users/bic/Desktop/cov-19_data/test_code/clean_test2h.txt'
    # target_amplicon = ['cov2_13_195bp', 'cov2_5_194bp', 'cov2_9_309bp', 'cov2_4_273bp', 'cov2_10_394bp',
    #                    'ACTB_263bp', 'influ_A_244bp', 'HAV_161bp', 'HKU1_151bp']
    # housekeeping = 'ACTB_263bp'
    # NTC = 'barcode61'

    outfile_header = '\t'.join(['\n\n#barcode', 'Result', 'read', 'base'] + target_amplicon) + '\n'
    # print(outfile_header)

    df = pd.read_table(updated_result)
    df = df.fillna('0/0/0')
    # print(df)

    summarized_result = [outfile_header]

    # get data of NTC
    ntc_record_dir = {}
    if NTC is not None:
        bc_df = df[df['#barcode'].str.match(NTC)]
        if bc_df.empty is True:
            logging.info("WARNING!   can not find NTC barcode data, ignore NTC in below result.")
            for amplicon in target_amplicon:
                ntc_record_dir[amplicon] = 0
        else:
            for amplicon in target_amplicon:
                part_df = bc_df[amplicon]
                part_df = part_df.str.split('/', expand=True).astype(float)
                # print(amplicon, "\n", part_df)
                ntc_record_dir[amplicon] = part_df[2].astype(int).values[0]
    else:
        for amplicon in target_amplicon:
            ntc_record_dir[amplicon] = 0

    barcode_list = df['#barcode'].unique()
    # barcode_list = ['barcode15']
    for bc in barcode_list:
        bc_df = df[df['#barcode'].str.match(bc)]
        # print(bc_df)
        record_dir = {}
        this_result_end = []

        for amplicon in target_amplicon:
            part_df = bc_df[amplicon]
            part_df = part_df.str.split('/', expand=True).astype(float)
            # print(amplicon, "\n", part_df)
            record_dir[amplicon] = part_df[2].astype(int).values[0]
            this_result_end.append(record_dir[amplicon])

        if housekeeping is None:
            housekeeping = "no_hp"
            record_dir[housekeeping] = 0

        cov2_50_plus = 0
        cov2_20_50 = 0
        cov2_10_20 = 0
        cov2_5_10 = 0
        cov2_1_5 = 0

        for key in record_dir:
            if key.startswith('cov2'):
                if int(record_dir[key] - ntc_record_dir[key]) >= 50:
                    cov2_50_plus += 1
                elif 20 <= int(record_dir[key] - ntc_record_dir[key]) < 50:
                    cov2_20_50 += 1
                elif 10 <= int(record_dir[key] - ntc_record_dir[key]) < 20:
                    cov2_10_20 += 1
                elif 5 <= int(record_dir[key] - ntc_record_dir[key]) < 10:
                    cov2_5_10 += 1
                elif 1 <= int(record_dir[key] - ntc_record_dir[key]) < 5:
                    cov2_1_5 += 1

        if cov2_50_plus >= 3:
            result_mark = 'POS_3'
        elif cov2_50_plus == 2:
            result_mark = 'POS_2'
        elif cov2_50_plus == 1 and int(cov2_20_50 + cov2_10_20 + cov2_5_10) >= 2:
            result_mark = 'POS_2'
        elif int(cov2_50_plus + cov2_20_50) >= 1:
            result_mark = 'POS_1'
        elif int(cov2_50_plus + cov2_20_50 + cov2_10_20 + cov2_5_10) >= 2:
            result_mark = 'POS_1'
        elif int(cov2_50_plus + cov2_20_50 + cov2_10_20 + cov2_5_10 + cov2_1_5) >= 3:
            result_mark = 'POS_1'
        elif int(cov2_50_plus + cov2_20_50 + cov2_10_20 + cov2_5_10 + cov2_1_5) >= 1:
            result_mark = 'POS_0'
        else:
            if int(record_dir[housekeeping] - ntc_record_dir[housekeeping]) >= 1000:
                result_mark = 'NEG'
            else:
                result_mark = 'UNK'
        # print(bc, '\t', result_mark)

        read_number = bc_df['read_number'].values[0]
        base_number = bc_df['total_base'].values[0]
        this_result = [bc, result_mark, read_number, base_number] + this_result_end
        this_result = '\t'.join([str(s) for s in this_result]) + '\n'
        summarized_result.append(this_result)

    with open(updated_result, 'a') as outfile:
        for line in summarized_result:
            outfile.writelines(line)
            logging.info(line.strip())
            # print(line.strip())


# summarize()
