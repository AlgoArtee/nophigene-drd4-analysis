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

  const preprocessForm = document.querySelector('[data-preprocess-form]');
  if (preprocessForm) {
    preprocessForm.addEventListener('submit', (event) => {
      const action = event.submitter?.value || '';
      if (!['find_region', 'select_methylation'].includes(action)) return;
      const status = document.querySelector('[data-preprocessing-status]');
      if (!status) return;
      const fill = status.querySelector('[data-preprocessing-progress-fill]');
      const track = fill?.parentElement;
      const caption = status.querySelector('[data-preprocessing-progress-caption]');
      const pill = status.querySelector('[data-preprocessing-status-pill]');
      const live = status.querySelector('[data-preprocessing-live]');
      const target = action === 'find_region' ? 48 : 94;
      const message = action === 'find_region'
        ? 'Resolving the gene interval…'
        : 'Filtering and saving the prepared manifest…';
      status.setAttribute('aria-busy', 'true');
      if (fill) fill.style.height = `${target}%`;
      if (track) track.setAttribute('aria-valuenow', String(target));
      if (caption) caption.textContent = 'Working…';
      if (pill) {
        pill.textContent = 'Working…';
        pill.classList.add('ready');
      }
      if (live) live.textContent = message;
    });
  }

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

  document.querySelectorAll('[data-literature-browser]').forEach((browser) => {
    const dataElement = browser.querySelector('[data-literature-data]');
    let findings = [];
    try {
      findings = JSON.parse(dataElement?.textContent || '[]');
    } catch (_error) {
      findings = [];
    }
    const search = browser.querySelector('[data-literature-search]');
    const priority = browser.querySelector('[data-literature-priority]');
    const source = browser.querySelector('[data-literature-source]');
    const variant = browser.querySelector('[data-literature-variant]');
    const pageSize = browser.querySelector('[data-literature-page-size]');
    const rows = browser.querySelector('[data-literature-rows]');
    const status = browser.querySelector('[data-literature-status]');
    const pageLabel = browser.querySelector('[data-literature-page]');
    const previous = browser.querySelector('[data-literature-previous]');
    const next = browser.querySelector('[data-literature-next]');
    let page = 1;

    const sources = [...new Set(findings.flatMap((item) => item.sources || [item.source_key]).filter(Boolean))].sort();
    const variants = [...new Set(findings.map((item) => item.variant).filter(Boolean))].sort();
    sources.forEach((value) => source.appendChild(new Option(value, value)));
    variants.forEach((value) => variant.appendChild(new Option(value, value)));

    const appendTextCell = (row, value) => {
      const cell = document.createElement('td');
      cell.textContent = value || '—';
      row.appendChild(cell);
      return cell;
    };
    const render = () => {
      const query = search.value.trim().toLocaleLowerCase();
      const filtered = findings.filter((item) => {
        if (priority.value && String(item.priority_tier) !== priority.value) return false;
        if (source.value && !(item.sources || [item.source_key]).includes(source.value)) return false;
        if (variant.value && item.variant !== variant.value) return false;
        if (!query) return true;
        return [item.finding, item.paper, item.title, item.phenotype, item.genotypes, item.variant,
          item.pmid, item.pmcid, item.doi, ...(item.sources || []), ...(item.evidence_tags || [])]
          .filter(Boolean).join(' ').toLocaleLowerCase().includes(query);
      });
      const size = Number(pageSize.value) || 20;
      const pages = Math.max(1, Math.ceil(filtered.length / size));
      page = Math.min(page, pages);
      const start = (page - 1) * size;
      const visible = filtered.slice(start, start + size);
      rows.replaceChildren();
      visible.forEach((item) => {
        const row = document.createElement('tr');
        appendTextCell(row, item.rank);
        const evidenceCell = appendTextCell(row, '');
        const priorityBadge = document.createElement('span');
        priorityBadge.className = `literature-badge priority-${item.priority_tier || 6}`;
        priorityBadge.textContent = item.priority_label || 'Gene-relevant evidence';
        evidenceCell.replaceChildren(priorityBadge);
        (item.evidence_tags || []).forEach((tag) => {
          const badge = document.createElement('span');
          badge.className = 'literature-tag';
          badge.textContent = tag.replaceAll('_', ' ');
          evidenceCell.appendChild(badge);
        });
        appendTextCell(row, item.finding || item.summary || 'Citation metadata only; no finding was synthesized.');
        const paperCell = appendTextCell(row, '');
        const paperText = item.paper || item.title || 'Publication';
        if (item.url) {
          const link = document.createElement('a');
          link.href = item.url;
          link.target = '_blank';
          link.rel = 'noopener noreferrer';
          link.textContent = paperText;
          paperCell.replaceChildren(link);
        } else {
          paperCell.textContent = paperText;
        }
        appendTextCell(row, item.phenotype);
        appendTextCell(row, item.genotypes);
        appendTextCell(row, item.variant);
        appendTextCell(row, (item.sources || [item.source_key]).filter(Boolean).join(', '));
        appendTextCell(row, [item.pmid && `PMID ${item.pmid}`, item.pmcid, item.doi && `DOI ${item.doi}`].filter(Boolean).join(' · '));
        rows.appendChild(row);
      });
      if (!visible.length) {
        const row = document.createElement('tr');
        const cell = document.createElement('td');
        cell.colSpan = 9;
        cell.className = 'empty';
        cell.textContent = 'No findings match these filters.';
        row.appendChild(cell);
        rows.appendChild(row);
      }
      status.textContent = filtered.length === findings.length
        ? `${findings.length} ranked finding${findings.length === 1 ? '' : 's'}`
        : `${filtered.length} of ${findings.length} findings match`;
      pageLabel.textContent = `Page ${page} of ${pages}`;
      previous.disabled = page <= 1;
      next.disabled = page >= pages;
    };
    [search, priority, source, variant, pageSize].forEach((control) => control.addEventListener('input', () => {
      page = 1;
      render();
    }));
    previous.addEventListener('click', () => {
      page = Math.max(1, page - 1);
      render();
    });
    next.addEventListener('click', () => {
      page += 1;
      render();
    });
    render();
  });

  const refreshPersonalStatistics = async () => {
    const shells = document.querySelectorAll('[data-personal-statistics]');
    if (!shells.length) return;
    try {
      const payload = await apiJson('/api/v2/personal-statistics');
      shells.forEach((shell) => {
        shell.querySelector('[data-personal-statistics-status]').textContent = `${payload.gene_count} gene${payload.gene_count === 1 ? '' : 's'}`;
        shell.querySelector('[data-personal-gene-count]').textContent = payload.gene_count;
        shell.querySelector('[data-personal-variant-count]').textContent = payload.totals.unique_variant_loci;
        shell.querySelector('[data-personal-probe-count]').textContent = payload.totals.unique_methylation_probes;
        const rows = shell.querySelector('[data-personal-statistics-rows]');
        rows.replaceChildren();
        payload.genes.forEach((item) => {
          const row = document.createElement('tr');
          [item.gene, item.variant_count, item.promoter_variant_count, item.gene_body_variant_count,
            item.methylation_probe_count, numberText(item.mean_beta), numberText(item.median_beta),
            item.genome_build, item.analysis_scope?.replaceAll('_', ' '), item.run_id].forEach((value) => addCell(row, value));
          rows.appendChild(row);
        });
        shell.querySelector('[data-personal-statistics-empty]').hidden = payload.gene_count !== 0;
        shell.querySelector('.table-shell').hidden = payload.gene_count === 0;
      });
    } catch (error) {
      shells.forEach((shell) => {
        const status = shell.querySelector('[data-personal-statistics-status]');
        status.textContent = 'Unavailable';
        status.title = error.message;
      });
    }
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

  const modelState = {settings: null, metadata: null};
  const modelStatusLabel = (value) => String(value || 'unavailable').replaceAll('_', ' ');

  const refreshModelSettings = async () => {
    const payload = await apiJson('/api/v2/model-settings');
    modelState.settings = payload;
    const credentialStatus = payload.alphagenome?.credential?.status || 'missing';
    const workerStatus = payload.alphagenome?.worker?.worker?.status || 'unavailable';
    document.querySelectorAll('[data-alphagenome-credential-status]').forEach((item) => {
      item.textContent = modelStatusLabel(credentialStatus);
      item.classList.toggle('ready', credentialStatus === 'verified');
    });
    document.querySelectorAll('[data-alphagenome-readiness]').forEach((item) => {
      const ready = credentialStatus === 'verified' && workerStatus === 'ready';
      item.textContent = ready ? 'Ready' : `${modelStatusLabel(credentialStatus)} credential · ${modelStatusLabel(workerStatus)} worker`;
      item.classList.toggle('ready', ready);
    });
    document.querySelectorAll('[data-model-catalog]').forEach((container) => {
      container.replaceChildren();
      payload.models.forEach((model) => {
        const row = document.createElement('article');
        row.className = 'model-catalog-row';
        const title = document.createElement('strong');
        title.textContent = model.name || model.id;
        const status = document.createElement('span');
        status.className = 'evidence-badge';
        let availability = model.id === 'alphagenome-api'
          ? (credentialStatus === 'verified' && workerStatus === 'ready' ? 'ready' : 'configuration required')
          : (['unsupported_for_epic_5mc', 'blocked_pending_verified_release', 'framework_only'].includes(model.status)
            ? 'scientifically blocked' : 'missing adapter/assets');
        status.textContent = availability;
        const detail = document.createElement('small');
        detail.textContent = model.description || model.scientific_gate || model.status || '';
        row.append(title, status, detail);
        container.appendChild(row);
      });
    });
    return payload;
  };

  const pollModelJob = async (jobId, onUpdate) => {
    for (;;) {
      const job = await apiJson(`/api/v2/model-jobs/${encodeURIComponent(jobId)}`);
      onUpdate?.(job);
      if (['succeeded', 'partial', 'failed', 'blocked', 'cancelled'].includes(job.status)) return job;
      await new Promise((resolve) => setTimeout(resolve, 2000));
    }
  };

  document.querySelectorAll('[data-alphagenome-credential-settings]').forEach((settings) => {
    const input = settings.querySelector('[data-alphagenome-api-key]');
    const message = settings.querySelector('[data-alphagenome-credential-message]');
    settings.querySelector('[data-alphagenome-credential-save]')?.addEventListener('click', async () => {
      const apiKey = input.value.trim();
      if (!apiKey) return setMessage(message, 'Enter an API key before saving.', 'error');
      try {
        await apiJson('/api/v2/model-settings/alphagenome-api/credential', {
          method: 'PUT', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({api_key: apiKey}),
        });
        input.value = '';
        setMessage(message, 'Credential encrypted and stored. Verify metadata access before running the model.', 'success');
        await refreshModelSettings();
      } catch (error) {
        setMessage(message, error.message, 'error');
      }
    });
    settings.querySelector('[data-alphagenome-credential-verify]')?.addEventListener('click', async () => {
      try {
        const job = await apiJson('/api/v2/model-settings/alphagenome-api/credential/verify', {method: 'POST', headers: {'Content-Type': 'application/json'}, body: '{}'});
        setMessage(message, `Verification queued (${job.id}). No personal data is sent.`, 'success');
        const finalState = await pollModelJob(job.id, (state) => setMessage(message, `Metadata verification: ${modelStatusLabel(state.stage)} · ${state.progress_percent}%`, state.status === 'failed' ? 'error' : 'success'));
        await refreshModelSettings();
        setMessage(message, finalState.status === 'succeeded' ? 'Credential verified and ontology metadata loaded.' : `Verification ${modelStatusLabel(finalState.status)}.`, finalState.status === 'succeeded' ? 'success' : 'error');
      } catch (error) {
        setMessage(message, error.message, 'error');
      }
    });
    settings.querySelector('[data-alphagenome-credential-delete]')?.addEventListener('click', async () => {
      try {
        await apiJson('/api/v2/model-settings/alphagenome-api/credential', {method: 'DELETE'});
        input.value = '';
        setMessage(message, 'Credential and cached provider metadata removed.', 'success');
        await refreshModelSettings();
      } catch (error) {
        setMessage(message, error.message, 'error');
      }
    });
  });

  const predictionRequestSettings = (step) => ({
    gene: step.dataset.gene || '',
    ontology_terms: Array.from(step._selectedOntology || []),
    modalities: Array.from(step.querySelectorAll('[data-alphagenome-modality]:checked')).map((item) => item.value),
    sequence_length: Number(step.querySelector('[data-alphagenome-sequence-length]')?.value || 1048576),
  });

  document.querySelectorAll('[data-model-step]').forEach((step) => {
    const runId = step.dataset.runId;
    if (!runId) return;
    step._selectedOntology = new Set();
    let previewDigest = '';
    const search = step.querySelector('[data-alphagenome-tissue-search]');
    const results = step.querySelector('[data-alphagenome-tissue-results]');
    const selection = step.querySelector('[data-alphagenome-tissue-selection]');
    const previewMessage = step.querySelector('[data-prediction-preview-message]');
    const previewPanel = step.querySelector('[data-prediction-preview-panel]');
    const consent = step.querySelector('[data-prediction-consent]');
    const submit = step.querySelector('[data-prediction-submit]');
    const jobMessage = step.querySelector('[data-prediction-job-message]');
    const renderSelection = () => {
      selection.replaceChildren();
      step._selectedOntology.forEach((curie) => {
        const button = document.createElement('button');
        button.type = 'button';
        button.className = 'ontology-chip';
        button.textContent = `${curie} ×`;
        button.addEventListener('click', () => { step._selectedOntology.delete(curie); renderSelection(); previewDigest = ''; submit.disabled = true; });
        selection.appendChild(button);
      });
    };
    let searchTimer;
    search?.addEventListener('input', () => {
      clearTimeout(searchTimer);
      searchTimer = setTimeout(async () => {
        const query = search.value.trim();
        results.replaceChildren();
        if (query.length < 2) return;
        try {
          const payload = await apiJson(`/api/v2/models/alphagenome-api/metadata?q=${encodeURIComponent(query)}`);
          payload.ontology_terms.slice(0, 20).forEach((term) => {
            const button = document.createElement('button');
            button.type = 'button';
            button.className = 'ontology-result';
            button.textContent = `${term.biosample_name || term.ontology_curie} · ${term.ontology_curie}`;
            button.disabled = step._selectedOntology.has(term.ontology_curie) || step._selectedOntology.size >= 5;
            button.addEventListener('click', () => {
              step._selectedOntology.add(term.ontology_curie);
              search.value = '';
              results.replaceChildren();
              renderSelection();
              previewDigest = '';
              submit.disabled = true;
            });
            results.appendChild(button);
          });
          if (!payload.ontology_terms.length) results.textContent = payload.reason || 'No matching ontology-backed biosamples.';
        } catch (error) {
          results.textContent = error.message;
        }
      }, 250);
    });
    step.querySelectorAll('[data-alphagenome-modality], [data-alphagenome-sequence-length]').forEach((control) => control.addEventListener('change', () => { previewDigest = ''; submit.disabled = true; }));
    step.querySelector('[data-prediction-preview]')?.addEventListener('click', async () => {
      previewPanel.hidden = true;
      try {
        const payload = await apiJson(`/api/v2/runs/${encodeURIComponent(runId)}/predictions/preview`, {
          method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(predictionRequestSettings(step)),
        });
        previewPanel.hidden = false;
        previewDigest = payload.payload_sha256;
        step.querySelector('[data-preview-selected]').textContent = payload.selected_variant_count;
        step.querySelector('[data-preview-omitted]').textContent = payload.omitted_variant_count;
        step.querySelector('[data-preview-terms]').textContent = payload.payload.ontology_terms.length;
        const rows = step.querySelector('[data-prediction-preview-rows]');
        rows.replaceChildren();
        payload.payload.variants.forEach((variant) => {
          const row = document.createElement('tr');
          const interval = variant.model_interval_0_based_half_open;
          [variant.variant, `${interval.start}–${interval.end} (${interval.length} bp)`, variant.reference_allele_verified ? 'yes' : 'no'].forEach((value) => addCell(row, value));
          rows.appendChild(row);
        });
        step.querySelector('[data-prediction-transfer-notice]').textContent = payload.external_transfer_notice;
        step.querySelector('[data-prediction-digest]').textContent = `SHA-256 ${payload.payload_sha256}`;
        consent.checked = false;
        submit.disabled = true;
        const blockers = payload.blockers || [];
        setMessage(previewMessage, payload.status === 'ready' ? 'Preflight passed. Review the complete payload below.' : `Blocked: ${blockers.join(', ')}`, payload.status === 'ready' ? 'success' : 'error');
      } catch (error) {
        previewDigest = '';
        setMessage(previewMessage, error.message, 'error');
      }
    });
    consent?.addEventListener('change', () => { submit.disabled = !(consent.checked && previewDigest); });
    submit?.addEventListener('click', async () => {
      try {
        const job = await apiJson(`/api/v2/runs/${encodeURIComponent(runId)}/predictions`, {
          method: 'POST', headers: {'Content-Type': 'application/json'},
          body: JSON.stringify({...predictionRequestSettings(step), payload_sha256: previewDigest, external_transfer_consent: true}),
        });
        submit.disabled = true;
        setMessage(jobMessage, `Queued model job ${job.id}.`, 'success');
        const finalState = await pollModelJob(job.id, (state) => setMessage(jobMessage, `${modelStatusLabel(state.stage)} · ${state.progress_percent}%`, state.status === 'failed' ? 'error' : 'success'));
        setMessage(jobMessage, `Model job ${modelStatusLabel(finalState.status)}. Open Predictions to inspect every available score or failure.`, ['succeeded', 'partial'].includes(finalState.status) ? 'success' : 'error');
        document.querySelectorAll(`[data-predictions-result][data-run-id="${CSS.escape(runId)}"]`).forEach((panel) => panel._refreshPredictions?.());
      } catch (error) {
        setMessage(jobMessage, error.message, 'error');
      }
    });
  });

  document.querySelectorAll('[data-predictions-result]').forEach((panel) => {
    const runId = panel.dataset.runId;
    const gene = panel.dataset.gene;
    let allScores = [];
    let page = 1;
    const filter = panel.querySelector('[data-prediction-score-filter]');
    const pageSize = panel.querySelector('[data-prediction-score-page-size]');
    const renderScores = () => {
      const query = (filter.value || '').trim().toLowerCase();
      const visible = allScores.filter((score) => !query || JSON.stringify(score).toLowerCase().includes(query));
      const size = Number(pageSize.value || 20);
      const pages = Math.max(1, Math.ceil(visible.length / size));
      page = Math.min(page, pages);
      const rows = panel.querySelector('[data-prediction-score-rows]');
      rows.replaceChildren();
      visible.slice((page - 1) * size, page * size).forEach((score) => {
        const row = document.createElement('tr');
        [score.variant, score.output_type || score.output_name, numberText(score.raw_score), numberText(score.quantile_score), score.scorer, score.track_name, score.gene_name || score.gene_id, [score.ontology_curie, score.biosample_name].filter(Boolean).join(' · ')].forEach((value) => addCell(row, value));
        rows.appendChild(row);
      });
      if (!rows.children.length) { const row = document.createElement('tr'); const cell = document.createElement('td'); cell.colSpan = 8; cell.textContent = 'No model scores match the current filter.'; row.appendChild(cell); rows.appendChild(row); }
      panel.querySelector('[data-prediction-score-page]').textContent = `Page ${page} of ${pages} · ${visible.length} rows`;
      panel.querySelector('[data-prediction-score-prev]').disabled = page <= 1;
      panel.querySelector('[data-prediction-score-next]').disabled = page >= pages;
    };
    filter.addEventListener('input', () => { page = 1; renderScores(); });
    pageSize.addEventListener('change', () => { page = 1; renderScores(); });
    panel.querySelector('[data-prediction-score-prev]').addEventListener('click', () => { page = Math.max(1, page - 1); renderScores(); });
    panel.querySelector('[data-prediction-score-next]').addEventListener('click', () => { page += 1; renderScores(); });
    const refresh = async () => {
      const payload = await apiJson(`/api/v2/runs/${encodeURIComponent(runId)}/predictions?gene=${encodeURIComponent(gene)}`);
      panel.querySelector('[data-predictions-status]').textContent = modelStatusLabel(payload.status);
      panel.querySelector('[data-predictions-observed]').textContent = payload.counts?.observed_allele_count ?? 0;
      panel.querySelector('[data-predictions-native]').textContent = payload.counts?.source_native_annotation_count ?? 0;
      panel.querySelector('[data-predictions-runs]').textContent = payload.counts?.model_run_count ?? 0;
      panel.querySelector('[data-predictions-scores]').textContent = payload.counts?.completed_prediction_count ?? 0;
      const jobs = panel.querySelector('[data-predictions-job-rows]');
      jobs.replaceChildren();
      (payload.model_runs || []).forEach((job) => { const row = document.createElement('tr'); [job.id, job.model_id, modelStatusLabel(job.status), `${job.progress_percent || 0}%`, job.started_at, job.finished_at].forEach((value) => addCell(row, value)); jobs.appendChild(row); });
      if (!jobs.children.length) { const row = document.createElement('tr'); const cell = document.createElement('td'); cell.colSpan = 6; cell.textContent = 'No model jobs have been submitted for this run.'; row.appendChild(cell); jobs.appendChild(row); }
      allScores = payload.predictions || [];
      renderScores();
      const failures = panel.querySelector('[data-predictions-failures]');
      if ((payload.failures || []).length) setMessage(failures, `${payload.failures.length} model output failure(s): ${payload.failures.map((item) => item.code || item.message).join(', ')}`, 'error');
      else failures.hidden = true;
      const notices = panel.querySelector('[data-predictions-notices]');
      if ((payload.notices || []).length) setMessage(notices, payload.notices.map((item) => item.message).join(' '), 'warning');
      else notices.hidden = true;
      return payload;
    };
    panel._refreshPredictions = refresh;
    refresh().then((payload) => {
      if (['queued', 'running'].includes(payload.status)) {
        const timer = setInterval(async () => { try { const next = await refresh(); if (!['queued', 'running'].includes(next.status)) clearInterval(timer); } catch (_error) { clearInterval(timer); } }, 2500);
      }
    }).catch(() => {});
  });

  Promise.allSettled([refreshPersonalStatistics(), refreshDandelionDatasets(), refreshDandelionHistory(), refreshDandelionHealth(), refreshModelSettings()]);
  showTask(document.body.dataset.initialTask || 'run');
})();
