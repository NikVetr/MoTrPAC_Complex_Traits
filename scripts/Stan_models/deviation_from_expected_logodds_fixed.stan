data {
    int<lower=1> row_n;
    int<lower=1> col_n;
    int<lower=1> n;
    int<lower=1> colcat_n;
    int<lower=0> total[n];
    int<lower=1,upper=row_n> row_index[n];
    int<lower=1,upper=col_n> col_index[n];
    int<lower=1,upper=colcat_n> colcat_index[col_n];
    int<lower=0> row_count[n];
    int<lower=0> col_count[n];
    int<lower=0> cell_count[n];
    vector<lower=0, upper=1>[n] overall_prob;
}
transformed data {
    vector[n] overall_logodds = logit(overall_prob);
}
parameters {
    //biases in deviations terms
    real overall_bias;
    vector[row_n] raw_row_bias;
    vector[col_n] raw_col_bias;
    vector[colcat_n] raw_colcat_bias;
    vector[n] raw_cell_bias;
    real<lower=0> row_bias_sd;
    real<lower=0> col_bias_sd;
    real<lower=0> colcat_bias_sd;
    real<lower=0> cell_bias_sd;
}
transformed parameters {
    //incorporate bias
    vector[colcat_n] colcat_bias = raw_colcat_bias * colcat_bias_sd;
    vector[col_n] col_bias = raw_col_bias * col_bias_sd + colcat_bias[colcat_index];
    vector[row_n] row_bias = raw_row_bias * row_bias_sd;
    vector[n] cell_bias = raw_cell_bias * cell_bias_sd;

    //add to expected logodds w/ no bias
    vector[n] cell_logodds = overall_logodds + overall_bias + row_bias[row_index] + col_bias[col_index] + cell_bias;
}
model {
    //bias params
    overall_bias ~ std_normal();

    raw_colcat_bias ~ std_normal();
    colcat_bias_sd ~ std_normal();

    raw_col_bias ~ std_normal();
    col_bias_sd ~ std_normal();

    raw_row_bias ~ std_normal();
    row_bias_sd ~ std_normal();

    raw_cell_bias ~ std_normal();
    cell_bias_sd ~ std_normal();

    //likelihood
    cell_count ~ binomial_logit(row_count, cell_logodds);

}
generated quantities {
    vector[n] cell_total_prob_bias = inv_logit(cell_logodds) - overall_prob;
}