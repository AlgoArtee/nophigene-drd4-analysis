(() => {
  const showTask = (key) => {
    document.querySelectorAll('[data-task-target]').forEach((button) => {
      button.setAttribute('aria-selected', String(button.dataset.taskTarget === key));
    });
    document.querySelectorAll('[data-task-panel]').forEach((panel) => {
      const active = panel.dataset.taskPanel === key;
      panel.hidden = !active;
      panel.classList.toggle('active', active);
    });
  };
  document.querySelectorAll('[data-task-target]').forEach((button) => {
    button.addEventListener('click', () => showTask(button.dataset.taskTarget));
  });
  document.querySelectorAll('[data-result-target]').forEach((button) => {
    button.addEventListener('click', () => {
      const key = button.dataset.resultTarget;
      document.querySelectorAll('[data-result-target]').forEach((item) => {
        item.setAttribute('aria-selected', String(item.dataset.resultTarget === key));
      });
      document.querySelectorAll('[data-result-panel]').forEach((panel) => {
        const active = panel.dataset.resultPanel === key;
        panel.hidden = !active;
        panel.classList.toggle('active', active);
      });
    });
    button.addEventListener('keydown', (event) => {
      const tabs = Array.from(document.querySelectorAll('[data-result-target]'));
      const current = tabs.indexOf(button);
      let next = current;
      if (event.key === 'ArrowRight') next = (current + 1) % tabs.length;
      else if (event.key === 'ArrowLeft') next = (current - 1 + tabs.length) % tabs.length;
      else if (event.key === 'Home') next = 0;
      else if (event.key === 'End') next = tabs.length - 1;
      else return;
      event.preventDefault();
      tabs[next].focus();
      tabs[next].click();
    });
  });
  document.querySelectorAll('[data-interaction-explorer]').forEach((explorer) => {
    const runId = explorer.dataset.runId;
    const output = explorer.querySelector('[data-interaction-output]');
    explorer.querySelectorAll('[data-interaction-hop]').forEach((button) => {
      button.addEventListener('click', async () => {
        output.hidden = false;
        output.textContent = 'Loading source-backed interaction graph…';
        try {
          const response = await fetch(`/api/v2/runs/${encodeURIComponent(runId)}/interactions?hops=${button.dataset.interactionHop}&node_cap=150`);
          const payload = await response.json();
          if (!response.ok) throw new Error(payload.message || 'Interaction expansion failed');
          output.textContent = JSON.stringify(payload, null, 2);
        } catch (error) {
          output.textContent = error.message;
        }
      });
    });
  });

  const apiJson = async (url, options = {}) => {
    const response = await fetch(url, options);
    const payload = await response.json();
    if (!response.ok) throw new Error(payload.message || payload.error?.message || `Request failed (${response.status})`);
    return payload;
  };
  const setMessage = (element, text, kind = '') => {
    if (!element) return;
    element.hidden = false;
    element.textContent = text;
    element.className = `message ${kind}`.trim();
  };
  const addCell = (row, value) => {
    const cell = document.createElement('td');
    cell.textContent = value === null || value === undefined || value === '' ? '—' : String(value);
    row.appendChild(cell);
  };
  const numberText = (value) => {
    const number = Number(value);
    if (!Number.isFinite(number)) return '—';
    return number === 0 || Math.abs(number) >= 0.001 ? number.toPrecision(4) : number.toExponential(3);
  };

  const renderDandelionResult = (report) => {
    const shell = document.querySelector('[data-dandelion-result]');
    if (!shell) return;
    document.querySelector('[data-standard-result]').hidden = true;
    shell.hidden = false;
    shell.querySelector('[data-dandelion-result-title]').textContent = report.summary.title;
    shell.querySelector('[data-dandelion-result-context]').textContent = `${report.run.phenotype || 'phenotype not declared'} · ${report.run.assembly || 'build not declared'}`;
    const summary = shell.querySelector('[data-dandelion-summary]');
    summary.replaceChildren();
    [
      ['Significant candidates', report.summary.candidate_count],
      ['Pairs tested', report.summary.tested_count],
      ['Rows shown', report.summary.returned_record_count],
      ['Method', 'DANDELION 0.1.0'],
    ].forEach(([label, value]) => {
      const article = document.createElement('article');
      const span = document.createElement('span');
      const strong = document.createElement('strong');
      span.textContent = label;
      strong.textContent = value;
      article.append(span, strong);
      summary.appendChild(article);
    });
    const sections = report.sections;
    shell.querySelector('[data-dandelion-objective]').textContent = `${sections.objective_data.status.replaceAll('_', ' ')}: ${sections.objective_data.reason}`;
    shell.querySelector('[data-dandelion-literature]').textContent = `${sections.literature.status.replaceAll('_', ' ')}: ${sections.literature.reason}`;
    shell.querySelector('[data-dandelion-medical]').textContent = `${sections.medical.status.replaceAll('_', ' ')}: ${sections.medical.reason}`;
    shell.querySelector('[data-dandelion-predictions]').textContent = `${sections.predictions.status.replaceAll('_', ' ')}: ${sections.predictions.reason}`;
    shell.querySelector('[data-dandelion-statistics-note]').textContent = `${sections.statistics.test_family}. Showing at most 20 ranked rows; full normalized and RDS artifacts remain in Run Details.`;
    const statistics = shell.querySelector('[data-dandelion-statistics-rows]');
    statistics.replaceChildren();
    (sections.statistics.results || []).slice(0, 20).forEach((item) => {
      const row = document.createElement('tr');
      addCell(row, item.exposure);
      addCell(row, item.candidate_gene);
      addCell(row, numberText(item.trans_p_value));
      addCell(row, numberText(item.gene_association_p_value));
      addCell(row, numberText(item.p_value));
      addCell(row, numberText(item.q_value));
      addCell(row, item.significant ? 'yes' : 'no');
      statistics.appendChild(row);
    });
    shell.querySelector('[data-dandelion-interactions-note]').textContent = sections.interactions.score_direction;
    const interactions = shell.querySelector('[data-dandelion-interaction-rows]');
    interactions.replaceChildren();
    (sections.interactions.edges || []).slice(0, 20).forEach((item) => {
      const row = document.createElement('tr');
      addCell(row, item.source_gene);
      addCell(row, item.target_gene);
      addCell(row, item.edge_type);
      addCell(row, numberText(item.native_score));
      addCell(row, item.tissue || 'not declared');
      interactions.appendChild(row);
    });
    shell.querySelector('[data-dandelion-run-details]').textContent = JSON.stringify(report.run_details, null, 2);
    shell.querySelectorAll('[data-cohort-result-target]').forEach((button) => {
      button.setAttribute('aria-selected', String(button.dataset.cohortResultTarget === 'summary'));
    });
    shell.querySelectorAll('[data-cohort-result-panel]').forEach((panel) => {
      const active = panel.dataset.cohortResultPanel === 'summary';
      panel.hidden = !active;
      panel.classList.toggle('active', active);
    });
    showTask('results');
  };

  document.querySelectorAll('[data-cohort-result-target]').forEach((button) => {
    button.addEventListener('click', () => {
      const key = button.dataset.cohortResultTarget;
      document.querySelectorAll('[data-cohort-result-target]').forEach((item) => item.setAttribute('aria-selected', String(item === button)));
      document.querySelectorAll('[data-cohort-result-panel]').forEach((panel) => {
        const active = panel.dataset.cohortResultPanel === key;
        panel.hidden = !active;
        panel.classList.toggle('active', active);
      });
    });
    button.addEventListener('keydown', (event) => {
      const tabs = Array.from(document.querySelectorAll('[data-cohort-result-target]'));
      const current = tabs.indexOf(button);
      let next = current;
      if (event.key === 'ArrowRight') next = (current + 1) % tabs.length;
      else if (event.key === 'ArrowLeft') next = (current - 1 + tabs.length) % tabs.length;
      else if (event.key === 'Home') next = 0;
      else if (event.key === 'End') next = tabs.length - 1;
      else return;
      event.preventDefault();
      tabs[next].focus();
      tabs[next].click();
    });
  });

  const refreshDandelionDatasets = async () => {
    const payload = await apiJson('/api/v2/statistical-datasets');
    document.querySelectorAll('[data-dandelion-dataset-select]').forEach((select) => {
      select.replaceChildren(new Option(payload.count ? 'Select a registered dataset' : 'No dataset registered', ''));
      payload.datasets.forEach((dataset) => {
        select.appendChild(new Option(`${dataset.name} · ${dataset.exposure_type} · ${dataset.assembly}`, dataset.id));
      });
    });
    document.querySelectorAll('[data-dandelion-datasets]').forEach((container) => {
      container.replaceChildren();
      if (!payload.count) {
        const empty = document.createElement('p');
        empty.className = 'empty';
        empty.textContent = 'No cohort datasets have been registered.';
        container.appendChild(empty);
      }
      payload.datasets.forEach((dataset) => {
        const card = document.createElement('article');
        card.className = 'dataset-record';
        const title = document.createElement('strong');
        title.textContent = dataset.name;
        const details = document.createElement('small');
        details.textContent = `${dataset.phenotype} · ${dataset.exposure_type} · ${dataset.assembly} · ${dataset.storage_mode.replaceAll('_', ' ')}`;
        card.append(title, details);
        if (dataset.storage_mode === 'registered_path') {
          const copy = document.createElement('button');
          copy.type = 'button';
          copy.className = 'quiet';
          copy.textContent = 'Create encrypted managed copy';
          copy.addEventListener('click', async () => {
            copy.disabled = true;
            copy.textContent = 'Encrypting…';
            try {
              await apiJson(`/api/v2/statistical-datasets/${encodeURIComponent(dataset.id)}/managed-copy`, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: '{}'});
              await refreshDandelionDatasets();
            } catch (error) {
              copy.disabled = false;
              copy.textContent = error.message;
            }
          });
          card.appendChild(copy);
        }
        container.appendChild(card);
      });
    });
  };

  const refreshDandelionHistory = async () => {
    const payload = await apiJson('/api/v2/statistical-analyses');
    document.querySelectorAll('[data-dandelion-history]').forEach((container) => {
      container.replaceChildren();
      if (!payload.count) {
        const empty = document.createElement('p');
        empty.className = 'empty';
        empty.textContent = 'No cohort statistical analyses have been submitted.';
        container.appendChild(empty);
      }
      payload.analyses.forEach((analysis) => {
        const card = document.createElement('article');
        card.className = 'card';
        const title = document.createElement('h3');
        title.textContent = `DANDELION · ${analysis.status}`;
        const details = document.createElement('p');
        details.textContent = `${analysis.progress_percent}% · dataset ${analysis.dataset_id}`;
        card.append(title, details);
        if (analysis.status === 'completed') {
          const open = document.createElement('button');
          open.type = 'button';
          open.textContent = 'Open canonical result';
          open.addEventListener('click', async () => renderDandelionResult(await apiJson(analysis.result_url)));
          card.appendChild(open);
        }
        if (['queued', 'running'].includes(analysis.status)) {
          const cancel = document.createElement('button');
          cancel.type = 'button';
          cancel.className = 'quiet';
          cancel.textContent = 'Request cancellation';
          cancel.addEventListener('click', async () => {
            cancel.disabled = true;
            try {
              await apiJson(`/api/v2/statistical-analyses/${encodeURIComponent(analysis.id)}/cancel`, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: '{}'});
              cancel.textContent = 'Cancellation requested';
            } catch (error) {
              cancel.disabled = false;
              cancel.textContent = error.message;
            }
          });
          card.appendChild(cancel);
        }
        if (['failed', 'cancelled'].includes(analysis.status)) {
          const retry = document.createElement('button');
          retry.type = 'button';
          retry.className = 'secondary';
          retry.textContent = 'Retry as a new analysis';
          retry.addEventListener('click', async () => {
            retry.disabled = true;
            try {
              const submitted = await apiJson(`/api/v2/statistical-analyses/${encodeURIComponent(analysis.id)}/retry`, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: '{}'});
              await refreshDandelionHistory();
              pollDandelionAnalysis(submitted.id, null).catch(() => {});
            } catch (error) {
              retry.disabled = false;
              retry.textContent = error.message;
            }
          });
          card.appendChild(retry);
        }
        container.appendChild(card);
      });
    });
  };

  const refreshDandelionHealth = async () => {
    const payload = await apiJson('/api/v2/health');
    const method = payload.dandelion;
    const status = method?.worker?.status || 'unavailable';
    document.querySelectorAll('[data-dandelion-status]').forEach((item) => {
      item.textContent = status === 'ready' ? 'Offline worker ready' : `Worker ${status}`;
      item.classList.toggle('ready', status === 'ready');
    });
    document.querySelectorAll('[data-dandelion-health]').forEach((item) => {
      item.textContent = `Status: ${status} · queue: ${method?.queue_depth ?? 'unknown'} · network: ${method?.network_policy || 'offline'}`;
    });
  };

  const pollDandelionAnalysis = async (analysisId, message) => {
    for (;;) {
      const state = await apiJson(`/api/v2/statistical-analyses/${encodeURIComponent(analysisId)}`);
      setMessage(message, `${state.stage.replaceAll('_', ' ')} · ${state.progress_percent}%`, state.status === 'failed' ? 'error' : 'success');
      if (state.status === 'completed') {
        renderDandelionResult(await apiJson(state.result_url));
        await refreshDandelionHistory();
        return;
      }
      if (['failed', 'cancelled'].includes(state.status)) return;
      await new Promise((resolve) => setTimeout(resolve, 2000));
    }
  };

  document.querySelectorAll('[data-dandelion-dataset-form]').forEach((form) => {
    form.addEventListener('submit', async (event) => {
      event.preventDefault();
      const message = form.closest('.card').querySelector('[data-dandelion-dataset-message]');
      try {
        const manifest = JSON.parse(new FormData(form).get('manifest'));
        const dataset = await apiJson('/api/v2/statistical-datasets', {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(manifest)});
        setMessage(message, `Registered ${dataset.name} (${dataset.id}); ${dataset.applicability.status} context applicability.`, 'success');
        await refreshDandelionDatasets();
      } catch (error) {
        setMessage(message, error.message, 'error');
      }
    });
  });
  document.querySelectorAll('[data-dandelion-analysis-form]').forEach((form) => {
    form.addEventListener('submit', async (event) => {
      event.preventDefault();
      const values = Object.fromEntries(new FormData(form));
      const message = form.closest('.card').querySelector('[data-dandelion-analysis-message]');
      try {
        const analysis = await apiJson('/api/v2/statistical-analyses', {
          method: 'POST',
          headers: {'Content-Type': 'application/json'},
          body: JSON.stringify({method: 'dandelion', dataset_id: values.dataset_id, parameters: {
            target_fdr: Number(values.target_fdr), cis_window_bp: Number(values.cis_window_bp),
            gene_association_threshold: Number(values.gene_association_threshold), chunk_size: Number(values.chunk_size),
          }}),
        });
        setMessage(message, `Queued ${analysis.id}. Estimated peak memory: ${Math.round(analysis.resource_estimate.estimated_peak_bytes / 1048576)} MiB.`, 'success');
        pollDandelionAnalysis(analysis.id, message).catch((error) => setMessage(message, error.message, 'error'));
      } catch (error) {
        setMessage(message, error.message, 'error');
      }
    });
  });

  Promise.allSettled([refreshDandelionDatasets(), refreshDandelionHistory(), refreshDandelionHealth()]);
  showTask(document.body.dataset.initialTask || 'run');
})();
